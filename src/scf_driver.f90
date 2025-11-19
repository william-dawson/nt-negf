!> SCF Driver - Handles I/O and calls scf_module for computation
!> Reads matrices once, passes to module, writes results
program scf_driver
  !! SCF module
  use scf_module
  !! NTPoly includes
  use datatypesmodule
  use loadbalancermodule
  use matrixconversionmodule
  use permutationmodule
  use processgridmodule
  use psmatrixalgebramodule
  use psmatrixmodule
  !! MPI
  use mpi
  implicit none

  !! Command-line parameters
  character(len=256) :: data_dir
  integer :: max_iterations
  real(ntreal) :: tolerance

  !! File paths
  character(len=256) :: h_file, s_file, k_file, idx_file

  !! NTPoly matrices
  type(matrix_ps) :: H, S, K, K_out
  type(matrix_ps) :: H_reordered, S_reordered, K_reordered

  !! Ordering and region info
  integer, dimension(:), allocatable :: order
  type(permutation_t) :: perm, perm_inv
  type(region_info_t) :: regions
  real(ntreal) :: mu_dft

  !! SCF parameters and results
  type(scf_params_t) :: scf_params
  real(ntreal) :: mu_final
  logical :: converged
  integer :: num_iterations

  !! MPI variables
  integer :: ierr, prov
  integer :: ii

  !! Initialize MPI and process grid
  call mpi_init_thread(mpi_thread_serialized, prov, ierr)
  call constructprocessgrid(mpi_comm_world)

  !! Parse command-line arguments
  call parse_arguments(data_dir, max_iterations, tolerance)

  !! Construct file paths
  h_file = trim(data_dir) // '/hamiltonian_sparse.mtx'
  s_file = trim(data_dir) // '/overlap_sparse.mtx'
  k_file = trim(data_dir) // '/density_kernel_sparse.mtx'
  idx_file = trim(data_dir) // '/order.txt'

  if (isroot()) then
    write(*,*) "==========================================================="
    write(*,*) "         SCF Chemical Potential Search Driver"
    write(*,*) "==========================================================="
    write(*,*) "Data directory:", trim(data_dir)
    write(*,*)
  end if

  !! ========================================================================
  !! STEP 1: Read matrices from disk (ONCE)
  !! ========================================================================
  if (isroot()) write(*,*) "Reading matrices from disk..."
  call constructmatrixfrommatrixmarket(H, h_file)
  call constructmatrixfrommatrixmarket(S, s_file)
  call constructmatrixfrommatrixmarket(K, k_file)

  !! ========================================================================
  !! STEP 2: Read ordering and region information (ONCE)
  !! ========================================================================
  if (isroot()) write(*,*) "Reading ordering information..."
  call get_order(idx_file, order, regions, mu_dft)

  if (isroot()) then
    write(*,'(A,F12.8,A)') " DFT chemical potential: ", mu_dft, " Ha"
    write(*,'(A,5I6)') " Regions (mid,lef1,lef2,rig1,rig2): ", &
         regions%mid, regions%lef1, regions%lef2, regions%rig1, regions%rig2
    write(*,*)
  end if

  !! ========================================================================
  !! STEP 3: Apply permutation to reorder matrices (ONCE)
  !! ========================================================================
  if (isroot()) write(*,*) "Applying permutation to matrices..."
  call constructdefaultpermutation(perm, h%actual_matrix_dimension)
  call constructdefaultpermutation(perm_inv, h%actual_matrix_dimension)

  perm%index_lookup = order + 1
  do ii = 1, h%actual_matrix_dimension
    perm%reverse_index_lookup(ii) = perm%index_lookup(ii)
    perm_inv%index_lookup(perm%index_lookup(ii)) = ii
    perm_inv%reverse_index_lookup(ii) = perm_inv%index_lookup(ii)
  end do

  call permutematrix(H, H_reordered, perm)
  call permutematrix(S, S_reordered, perm)
  call permutematrix(K, K_reordered, perm)

  !! Snap to sparsity pattern
  call snapmatrixtosparsitypattern(H_reordered, K_reordered)
  call snapmatrixtosparsitypattern(S_reordered, K_reordered)

  !! ========================================================================
  !! STEP 4: Setup SCF parameters
  !! ========================================================================
  scf_params%max_iterations = max_iterations
  scf_params%tolerance = tolerance
  scf_params%bracket_width = 0.05d0
  scf_params%mu_initial = mu_dft

  if (isroot()) then
    write(*,*) "SCF Parameters:"
    write(*,'(A,I4)')      "  Max iterations: ", scf_params%max_iterations
    write(*,'(A,ES10.3)')  "  Tolerance:      ", scf_params%tolerance
    write(*,'(A,F10.6,A)') "  Bracket width:  +/- ", scf_params%bracket_width, " Ha"
  end if

  !! ========================================================================
  !! STEP 5: Call SCF solver (matrices stay in memory, no repeated I/O!)
  !! ========================================================================
  call scf_solve(H_reordered, S_reordered, K_reordered, regions, &
                 scf_params, mu_final, converged, num_iterations)

  !! ========================================================================
  !! STEP 6: Undo permutation and write results
  !! ========================================================================
  if (isroot()) then
    write(*,*)
    write(*,*) "==========================================================="
    write(*,*) "Writing results..."
  end if

  call permutematrix(K_reordered, K_out, perm_inv)
  call writematrixtomatrixmarket(K_out, "density_scf.mtx")

  !! Save SCF summary
  if (isroot()) then
    call save_scf_summary(mu_dft, mu_final, converged, num_iterations)
  end if

  !! Cleanup
  deallocate(order)
  call destructprocessgrid
  call mpi_finalize(ierr)

  if (isroot()) then
    write(*,*)
    write(*,*) "Done!"
    write(*,*) "==========================================================="
  end if

contains

  !> Parse command-line arguments
  subroutine parse_arguments(data_dir, max_iter, tol)
    character(len=*), intent(out) :: data_dir
    integer, intent(out) :: max_iter
    real(ntreal), intent(out) :: tol
    character(len=256) :: arg
    integer :: nargs, i

    !! Default values
    data_dir = './example/data-sw'
    max_iter = 15
    tol = 1.0d-3

    !! Parse arguments
    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
      call get_command_argument(i, arg)
      select case (trim(arg))
      case ('--data-dir')
        if (i+1 <= nargs) then
          i = i + 1
          call get_command_argument(i, data_dir)
        end if
      case ('--max-iter')
        if (i+1 <= nargs) then
          i = i + 1
          call get_command_argument(i, arg)
          read(arg, *) max_iter
        end if
      case ('--tolerance')
        if (i+1 <= nargs) then
          i = i + 1
          call get_command_argument(i, arg)
          read(arg, *) tol
        end if
      end select
      i = i + 1
    end do
  end subroutine parse_arguments

  !> Read ordering information from file
  subroutine get_order(idx_file, order, regions, mu)
    character(len=*), intent(in) :: idx_file
    integer, dimension(:), allocatable, intent(out) :: order
    type(region_info_t), intent(out) :: regions
    real(ntreal), intent(out) :: mu
    integer, parameter :: fh = 16
    integer :: length

    open(fh, file = idx_file, status = 'old')
    read(fh, *) length
    allocate(order(length))
    read(fh, *) order
    read(fh, *) regions%mid, regions%lef1, regions%lef2, regions%rig1, regions%rig2
    read(fh, *) mu
    close(fh)
  end subroutine get_order

  !> Save SCF summary to file
  subroutine save_scf_summary(mu_initial, mu_final, converged, niter)
    real(ntreal), intent(in) :: mu_initial, mu_final
    logical, intent(in) :: converged
    integer, intent(in) :: niter
    integer, parameter :: fh = 20

    open(fh, file='scf_summary.txt', status='replace')
    write(fh, '(A)') '# SCF Summary'
    write(fh, '(A,F12.8)') 'mu_initial_Ha = ', mu_initial
    write(fh, '(A,F12.8)') 'mu_final_Ha = ', mu_final
    write(fh, '(A,F12.8)') 'mu_shift_Ha = ', mu_final - mu_initial
    write(fh, '(A,L1)')    'converged = ', converged
    write(fh, '(A,I4)')    'iterations = ', niter
    close(fh)

    write(*,*)
    write(*,*) "SCF Summary:"
    write(*,'(A,F12.8,A)') "  Initial mu: ", mu_initial, " Ha"
    write(*,'(A,F12.8,A)') "  Final mu:   ", mu_final, " Ha"
    write(*,'(A,F12.8,A)') "  Shift:      ", mu_final - mu_initial, " Ha"
    write(*,'(A,L1)')      "  Converged:  ", converged
    write(*,'(A,I4)')      "  Iterations: ", niter
    write(*,*)
    write(*,*) "Results written to:"
    write(*,*) "  - density_scf.mtx"
    write(*,*) "  - scf_summary.txt"
  end subroutine save_scf_summary

end program scf_driver
