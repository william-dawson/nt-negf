!> Self-consistent chemical potential search using linear interpolation
!> Uses regula falsi (false position) method with bisection fallback
!> This is a pure Fortran implementation of the Python script scan_mu_scf.py
program scf_driver
  !! includes related to libNEGF
  use libnegf
  use lib_param
  !! includes related to NTPoly
  use datatypesmodule
  use loadbalancermodule
  use matrixconversionmodule
  use permutationmodule
  use processgridmodule
  use psmatrixalgebramodule
  use psmatrixmodule
  use smatrixalgebramodule
  use smatrixmodule
  use tripletlistmodule
  use tripletmodule
  !! any other modules
  use mpi
  implicit none

  !! SCF parameters
  integer :: max_iterations
  real(ntreal) :: tolerance
  real(ntreal) :: mu_low, mu_high, mu, mu_dft
  real(ntreal) :: bracket_width
  character(len=256) :: data_dir

  !! File paths
  character(len=256) :: h_file, s_file, k_file, idx_file

  !! Convergence tracking
  real(ntreal) :: target_trace, current_trace, trace_error
  real(ntreal) :: trace_low, trace_high  !! Traces at bracket endpoints
  real(ntreal), dimension(:), allocatable :: mu_history, trace_history
  integer :: iteration, num_iterations
  logical :: converged, have_trace_low, have_trace_high

  !! MPI variables
  integer :: ierr, prov

  !! Initialize MPI
  call mpi_init_thread(mpi_thread_serialized, prov, ierr)
  call constructprocessgrid(mpi_comm_world)

  !! Parse command-line arguments
  call parse_arguments(data_dir, max_iterations, tolerance)

  !! Construct file paths
  h_file = trim(data_dir) // '/hamiltonian_sparse.mtx'
  s_file = trim(data_dir) // '/overlap_sparse.mtx'
  k_file = trim(data_dir) // '/density_kernel_sparse.mtx'
  idx_file = trim(data_dir) // '/order.txt'

  !! Read DFT chemical potential from order.txt
  call read_mu_dft(idx_file, mu_dft)

  !! Print header
  if (isroot()) then
    write(*,*) "==========================================================="
    write(*,*) "SCF Chemical Potential Search (Linear Interpolation)"
    write(*,*) "==========================================================="
    write(*,*) "Data directory:", trim(data_dir)
    write(*,*) "Max iterations:", max_iterations
    write(*,*) "Tolerance:", tolerance
    write(*,'(A,F12.8,A)') " DFT chemical potential: ", mu_dft, " Ha"
    write(*,*)
  end if

  !! Set initial bracket centered on DFT mu
  bracket_width = 0.05d0
  mu_low = mu_dft - bracket_width
  mu_high = mu_dft + bracket_width

  if (isroot()) then
    write(*,'(A,F10.6,A,F10.6,A)') " Initial bracket: [", mu_low, ", ", mu_high, "] Ha"
    write(*,'(A,F10.6,A)') " Bracket width: +/- ", bracket_width, " Ha"
    write(*,*)
  end if

  !! Allocate history arrays
  allocate(mu_history(max_iterations))
  allocate(trace_history(max_iterations))

  !! Initialize target trace (will be set on first iteration)
  target_trace = 0.0d0
  converged = .false.
  num_iterations = 0
  have_trace_low = .false.
  have_trace_high = .false.

  !! SCF loop
  do iteration = 1, max_iterations
    num_iterations = iteration

    !! Choose next mu value based on linear interpolation or bisection
    if (have_trace_low .and. have_trace_high) then
      !! Linear interpolation (regula falsi / false position method)
      !! Since trace increases with mu, interpolate to find mu for target_trace
      if (abs(trace_high - trace_low) > 1.0d-12) then
        mu = mu_low + (mu_high - mu_low) * (target_trace - trace_low) / (trace_high - trace_low)
        !! Safety check: ensure mu is within bracket
        mu = max(mu_low, min(mu_high, mu))
      else
        !! Fall back to bisection if traces are too close
        mu = (mu_low + mu_high) / 2.0d0
      end if
    else
      !! First iterations: use bisection until we have both bracket traces
      mu = (mu_low + mu_high) / 2.0d0
    end if

    if (isroot()) then
      write(*,*) "=== Iteration", iteration, "==="
      write(*,'(A,F12.8,A)') " mu = ", mu, " Ha"
      if (have_trace_low .and. have_trace_high) then
        write(*,'(A,F10.6,A,F10.6,A)') " Bracket: [", mu_low, ", ", mu_high, "]"
        write(*,'(A,F10.4,A,F10.4,A)') " Traces:  [", trace_low, ", ", trace_high, "]"
      end if
    end if

    !! Run NEGF iteration with current mu
    call run_negf_iteration(h_file, s_file, k_file, idx_file, mu, &
                            target_trace, current_trace)

    !! Store history
    mu_history(iteration) = mu
    trace_history(iteration) = current_trace

    !! Compute error
    trace_error = abs(current_trace - target_trace)

    if (isroot()) then
      write(*,'(A,F12.6,A,F12.6,A,F12.6,A)') " Trace: ", current_trace, &
           " (target: ", target_trace, ", error: ", trace_error, ")"
    end if

    !! Check convergence
    if (trace_error < tolerance) then
      converged = .true.
      if (isroot()) then
        write(*,*)
        write(*,*) "*** CONVERGED in", iteration, "iterations ***"
        write(*,'(A,F12.8,A)') " Final mu = ", mu, " Ha"
      end if
      exit
    end if

    !! Update bracket (higher mu -> higher trace)
    if (current_trace < target_trace) then
      mu_low = mu
      trace_low = current_trace
      have_trace_low = .true.
      if (isroot()) write(*,*) "-> Too low, increase mu"
    else
      mu_high = mu
      trace_high = current_trace
      have_trace_high = .true.
      if (isroot()) write(*,*) "-> Too high, decrease mu"
    end if

    if (isroot()) write(*,*)
  end do

  !! Print final results
  if (.not. converged .and. isroot()) then
    write(*,*)
    write(*,*) "*** Max iterations reached ***"
    write(*,'(A,F12.8,A,F12.8)') " Final mu = ", mu, " Ha, error = ", trace_error
  end if

  !! Save results
  if (isroot()) then
    call save_scf_data(mu_history, trace_history, num_iterations)
  end if

  !! Cleanup
  deallocate(mu_history, trace_history)
  call destructprocessgrid
  call mpi_finalize(ierr)

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

  !> Run a single NEGF iteration with specified mu
  subroutine run_negf_iteration(h_file, s_file, k_file, idx_file, mu_val, &
                                 target_trace, current_trace)
    character(len=*), intent(in) :: h_file, s_file, k_file, idx_file
    real(ntreal), intent(in) :: mu_val
    real(ntreal), intent(inout) :: target_trace
    real(ntreal), intent(out) :: current_trace

    !! Local variables - matrices
    type(matrix_ps) :: H, S, K, Hc, Sc, Kc
    type(matrix_lsc) :: Hloc, Sloc, Kloc

    !! Local variables - libNEGF
    type(Tnegf), target :: pnegf
    type(lnParams) :: params
    real(kind(1.d0)), dimension(:,:), pointer :: transmission

    !! Local variables - ordering
    integer, dimension(:), allocatable :: order
    integer :: mid, lef1, lef2, rig1, rig2
    type(permutation_t) :: perm
    integer, dimension(2) :: surfstart, surfend, contend
    integer, dimension(1) :: plend
    integer, dimension(:), allocatable :: cblk
    real(ntreal) :: mu_dft

    !! Local variables - temporary
    integer :: ii, ierr
    type(tripletlist_c) :: tlist
    type(triplet_c) :: trip
    real(ntreal) :: dotv
    complex(ntcomplex) :: dotc

    !! Setup libNEGF
    call setup_negf_local(pnegf)

    !! Read in the BigDFT Matrices
    call constructmatrixfrommatrixmarket(H, h_file)
    call constructmatrixfrommatrixmarket(S, s_file)
    call constructmatrixfrommatrixmarket(K, k_file)

    !! Compute target trace from original K and S
    call dotmatrix(K, S, dotv)
    if (target_trace == 0.0d0) then
      target_trace = dotv
      if (isroot()) then
        write(*,'(A,F12.6)') " Target trace: ", target_trace
      end if
    end if

    !! Read in the indexing information
    call get_order(idx_file, order, mid, lef1, lef2, rig1, rig2, mu_dft)

    allocate(cblk(2))
    plend(1) = mid
    surfstart(1) = mid + 1
    surfstart(2) = lef2 + 1
    surfend(1) = mid
    surfend(2) = lef2
    contend(1) = lef2
    contend(2) = rig2
    cblk = (1, 1)

    !! Reorder the matrices
    call constructdefaultpermutation(perm, h%actual_matrix_dimension)
    perm%index_lookup = order + 1
    do ii = 1, h%actual_matrix_dimension
       perm%reverse_index_lookup(ii) = perm%index_lookup(ii)
    end do
    call permutematrix(H, H, perm)
    call permutematrix(S, S, perm)
    call permutematrix(K, K, perm)

    !! Snap to sparsity pattern
    call snapmatrixtosparsitypattern(H, K)
    call snapmatrixtosparsitypattern(S, K)

    !! Convert to complex
    call convertmatrixtocomplex(H, Hc)
    call convertmatrixtocomplex(S, Sc)
    call convertmatrixtocomplex(K, Kc)

    !! Gather locally
    call gathermatrixtoprocess(Hc, Hloc)
    call gathermatrixtoprocess(Sc, Sloc)
    call gathermatrixtoprocess(Kc, Kloc)

    !! Convert to libNEGF matrices
    call create_HS(pnegf, 1) ! 1 k-point
    call create_DM(pnegf, 1)
    call set_h(pnegf, Hloc%rows, Hloc%values, Hloc%inner_index, Hloc%outer_index + 1)
    call set_s(pnegf, Sloc%rows, Sloc%values, Sloc%inner_index, Sloc%outer_index + 1)

    !! Initialize libNEGF structure
    call init_contacts(pnegf, 2)
    call init_structure(pnegf, 2, surfstart, surfend, contend, 1, plend, cblk)

    !! Set libNEGF parameters with current mu
    call get_params(pnegf, params)
    params%mu(1:2) = mu_val
    params%Emin = mu_val - 2.0d0
    params%Emax = mu_val + 1.0d0
    params%Estep = 1.d-3
    params%kbT_dm(1:2) = 1e-3
    params%kbT_t(1:2) = 1e-3
    params%ec = -2.0
    call set_params(pnegf, params)

    !! Compute density
    call compute_density_dft(pnegf)
    call mpi_allreduce(mpi_in_place, pnegf%rho%nzval, &
         & size(pnegf%rho%nzval), MPI_DOUBLE_COMPLEX,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)

    !! Compute transmission
    call compute_current(pnegf)
    call associate_transmission(pnegf, transmission)
    call mpi_allreduce(mpi_in_place, transmission, &
         & size(transmission, 1) * size(transmission, 2), &
         & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    !! Write tunneling data (only on final iteration if desired)
    ! if (isroot()) call write_tunneling_and_dos(pnegf)

    !! Patch back in the libNEGF density update
    call constructtripletlist(tlist)
    if (isroot()) then
      call extract_lnegf(pnegf%rho%ncol, pnegf%rho%rowpnt,&
                         pnegf%rho%colind, pnegf%rho%nzval, &
                         mid, lef1, lef2, rig1, rig2, tlist, trip)

      call extract_ntpoly(kloc%columns, kloc%outer_index, &
                          kloc%inner_index, kloc%values, &
                          mid, lef1, lef2, rig1, rig2, tlist, trip)
    end if
    call constructemptymatrix(Kc, S)
    call fillmatrixfromtripletlist(Kc, tlist)
    call convertmatrixtoreal(Kc, K)

    !! Compute current trace
    call dotmatrix(K, S, current_trace)

    !! Cleanup
    deallocate(order, cblk)
    call destroy_negf(pnegf)

  end subroutine run_negf_iteration

  !> Read just the DFT chemical potential from order.txt
  subroutine read_mu_dft(idx_file, mu)
    character(len=*), intent(in) :: idx_file
    real(ntreal), intent(out) :: mu
    integer, parameter :: fh = 15
    character(len=1024) :: line

    open(fh, file = idx_file, status = 'old')
    read(fh, *) ! Skip line 1: length
    read(fh, *) ! Skip line 2: order array
    read(fh, *) ! Skip line 3: region boundaries
    read(fh, *) mu ! Line 4: chemical potential
    close(fh)
  end subroutine read_mu_dft

  !> Read in the file with the ordering.
  subroutine get_order(idx_file, order, mid, lef1, lef2, rig1, rig2, mu)
    character(len=*), intent(in) :: idx_file
    integer, dimension(:), allocatable, intent(out) :: order
    integer, intent(out) :: mid, lef1, lef2, rig1, rig2
    real(ntreal), intent(out) :: mu
    integer, parameter :: fh = 16
    integer :: ii, length

    open(fh, file = idx_file, status = 'old')
    read(fh, *) length
    allocate(order(length))
    read(fh, *) order
    read(fh, *) mid, lef1, lef2, rig1, rig2
    read(fh, *) mu
    close(fh)
  end subroutine get_order

  !> Test if a row / column is in the appropriate region.
  logical function in_region(col, row, mid, lef1, lef2, rig1, rig2)
    integer, intent(in) :: col, row, mid, lef1, lef2, rig1, rig2
    in_region = (col <= mid .and. row <= mid) .or. &
        (col > mid .and. col <= lef1 .and. row <= mid) .or. &
        (row > mid .and. row <= lef1 .and. col <= mid) .or. &
        (col > lef2 .and. col <= rig1 .and. row <= mid) .or. &
        (row > lef2 .and. row <= rig1 .and. col <= mid)
  end function in_region

  !> Extract elements from libNEGF density in specific regions
  subroutine extract_lnegf(ncol, rowpnt, colind, values, &
                           mid, lef1, lef2, rig1, rig2, tlist, trip)
    integer, intent(in) :: ncol
    integer, intent(in) :: rowpnt(:), colind(:)
    complex(kind=8), intent(in) :: values(:)
    integer, intent(in) :: mid, lef1, lef2, rig1, rig2
    type(tripletlist_c), intent(inout) :: tlist
    type(triplet_c), intent(inout) :: trip
    integer :: jj, kk, row

    do jj = 1, ncol
      do kk = rowpnt(jj), rowpnt(jj+1) - 1
        row = colind(kk)
        if (in_region(jj, row, mid, lef1, lef2, rig1, rig2)) then
          trip%index_column = jj
          trip%index_row    = row
          trip%point_value  = values(kk)
          call appendtotripletlist(tlist, trip)
        end if
      end do
    end do
  end subroutine extract_lnegf

  !> Extract elements from NTPoly density not in libNEGF regions
  subroutine extract_ntpoly(ncol, rowpnt, colind, values, &
                            mid, lef1, lef2, rig1, rig2, tlist, trip)
    integer, intent(in) :: ncol
    integer, intent(in) :: rowpnt(:), colind(:)
    complex(kind=8), intent(in) :: values(:)
    integer, intent(in) :: mid, lef1, lef2, rig1, rig2
    type(tripletlist_c), intent(inout) :: tlist
    type(triplet_c), intent(inout) :: trip
    integer :: jj, kk, row

    do jj = 1, ncol
      do kk = rowpnt(jj) + 1, rowpnt(jj+1)
        row = colind(kk)
        if (.not. in_region(jj, row, mid, lef1, lef2, rig1, rig2)) then
          trip%index_column = jj
          trip%index_row    = row
          trip%point_value  = values(kk)
          call appendtotripletlist(tlist, trip)
        end if
      end do
    end do
  end subroutine extract_ntpoly

  !> Save SCF iteration data to file
  subroutine save_scf_data(mu_hist, trace_hist, niter)
    real(ntreal), dimension(:), intent(in) :: mu_hist, trace_hist
    integer, intent(in) :: niter
    integer, parameter :: fh = 20
    integer :: i

    open(fh, file='scf_data.txt', status='replace')
    write(fh, '(A)') '# Iteration  mu(Ha)  trace'
    do i = 1, niter
      write(fh, '(I3,2X,F12.8,2X,F12.6)') i, mu_hist(i), trace_hist(i)
    end do
    close(fh)

    write(*,*)
    write(*,*) "Data saved to: scf_data.txt"
  end subroutine save_scf_data

  !> Setup libNEGF with proper MPI interface
  subroutine setup_negf_local(pnegf)
    use mpi_f08
    type(Tnegf) :: pnegf
    call init_negf(pnegf)
    call set_mpi_bare_comm(pnegf, MPI_COMM_WORLD)
  end subroutine setup_negf_local

end program scf_driver
