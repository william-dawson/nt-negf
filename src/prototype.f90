!> A prototype of in memory interaction between libNEGF and NTPoly
program prototype
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
  !! the matrices prepared by BigDFT
  type(matrix_ps) :: H, S, K, Hc, Sc, Kc
  type(matrix_lsc) :: Hloc, Sloc, Kloc
  character(len=80) :: h_file, s_file, k_file, idx_file, o_file
  real(ntreal) :: mu, mu_initial
  !! libnegf types
  type(Tnegf), target :: pnegf
  type(lnParams) :: params
  real(kind(1.d0)), dimension(:,:), pointer :: transmission
  !! ordering information
  integer, dimension(:), allocatable :: order
  integer :: mid, lef1, lef2, rig1, rig2
  type(permutation_t) :: perm
  integer, dimension(2) :: surfstart, surfend, contend
  integer, dimension(1) :: plend
  integer, dimension(:), allocatable :: cblk
  !! the update prepared by libNEGF
  type(matrix_ps) :: Kupdate
  !! temporary variables
  integer :: prov, ierr, ii, jj, kk, col, row
  character(len=80) :: temps
  type(tripletlist_c) :: tlist
  type(triplet_c) :: trip
  real(ntreal) :: dotv
  complex(ntcomplex) :: dotc
  !! trace tracking variables
  real(ntreal) :: target_trace, current_trace

  !! Setup MPI, NTPoly, libNEGF
  call mpi_init_thread(mpi_thread_serialized, prov, ierr)
  call constructprocessgrid(mpi_comm_world)
  call setup_negf(pnegf)

  !! Read in the command line arguments
  call get_command_argument(1, h_file)
  call get_command_argument(2, s_file)
  call get_command_argument(3, k_file)
  call get_command_argument(4, idx_file)
  ! call get_command_argument(5, o_file)

  !! Read in the BigDFT Matrices
  call constructmatrixfrommatrixmarket(H, h_file)
  call constructmatrixfrommatrixmarket(S, s_file)
  call constructmatrixfrommatrixmarket(K, k_file)

  !! The number of occupied orbitals read from the trace
  call dotmatrix(K, S, dotv)
  if (isroot()) then
      write(*,*) "Target trace should be:", dotv
  end if
  target_trace = dotv

  !! Read in the indexing information
  call get_order(idx_file, order, mid, lef1, lef2, rig1, rig2, mu)
  write(*,*) mid, lef1, lef2, rig1, rig2
  if (isroot()) then
    write(*,*) "DFT chemical potential:", mu, "Ha"
  end if
  !! Adjust mu for open NEGF system to conserve charge
  !! The chemical potential shifts by ~+0.01 Ha when contacts are added
  mu_initial = mu + 0.01d0
  mu = mu_initial
  if (isroot()) then
    write(*,*) "Initial mu for NEGF:", mu, "Ha"
  end if

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
     perm%reverse_index_lookup(ii) = perm%index_lookup(II)
  end do
  call permutematrix(H, H, perm)
  call permutematrix(S, S, perm)
  call permutematrix(K, K, perm)
  ! call printmatrix(k)

  !! snap to sparsity pattern
  call snapmatrixtosparsitypattern(H, K)
  call snapmatrixtosparsitypattern(S, K)

  !! convert to complex
  call convertmatrixtocomplex(H, Hc)
  call convertmatrixtocomplex(S, Sc)
  call convertmatrixtocomplex(K, Kc)
  call dotmatrix(Kc, Sc, dotv)
  if (isroot()) then
      write(*,*) "Trace after sparsity switch", dotv
  end if

  !! Gather locally
  call gathermatrixtoprocess(Hc, Hloc)
  call gathermatrixtoprocess(Sc, Sloc)
  call gathermatrixtoprocess(Kc, Kloc)
  call dotmatrix(Kloc, Sloc, dotc)
  if (isroot()) then
      write(*,*) "Trace after local", dotc
  end if

  if (isroot()) then
      write(*,*) "Setup Done"
  end if

  !! Convert to libNEGF matrices
  call create_HS(pnegf, 1) ! 1 k-point
  call create_DM(pnegf, 1)
  call set_h(pnegf, Hloc%rows, Hloc%values, Hloc%inner_index, Hloc%outer_index + 1)
  call set_s(pnegf, Sloc%rows, Sloc%values, Sloc%inner_index, Sloc%outer_index + 1)

  !! Initialize libNEGF structure
  call init_contacts(pnegf, 2)
  call init_structure(pnegf, 2, surfstart, surfend, contend, 1, plend, cblk)

  !! Set libNEGF parameters
  call get_params(pnegf, params)
  params%mu(1:2) = mu
  params%Emin = mu - 2.0d0
  params%Emax = mu + 1.0d0
  params%Estep = 1.d-3
  params%kbT_dm(1:2) = 1e-3
  params%kbT_t(1:2) = 1e-3
  params%ec = -2.0
  if (isroot()) then
    write(*,*) "Energy window: Emin =", params%Emin, "Emax =", params%Emax
    write(*,*) "Chemical potential mu =", params%mu(1)
  end if
  call set_params(pnegf, params)

  !! Compute density
  call compute_density_dft(pnegf)
  call mpi_allreduce(mpi_in_place, pnegf%rho%nzval, &
       & size(pnegf%rho%nzval), MPI_DOUBLE_COMPLEX,&
       & MPI_SUM, MPI_COMM_WORLD, ierr)
  if (isroot()) write(*,*) "Density Done"

  !! Compute transmission
  call compute_current(pnegf)
  call associate_transmission(pnegf, transmission)
  call mpi_allreduce(mpi_in_place, transmission, &
       & size(transmission, 1) * size(transmission, 2), &
       & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (isroot()) then
    call write_tunneling_and_dos(pnegf)
    write(*,*) "Transmission Done"
  end if

  !! Patch back in the libNEGF density update
  call constructtripletlist(tlist)
  if (isroot()) then
    write(*,*) "Extracting elements from libNEGF and BigDFT..."
    call extract_lnegf(pnegf%rho%ncol, pnegf%rho%rowpnt,&
                       pnegf%rho%colind, pnegf%rho%nzval)

    call extract_ntpoly(kloc%columns, kloc%outer_index, &
                        kloc%inner_index, kloc%values)
    write(*,*) "Total triplets in combined list:", tlist%currentsize
  end if
  call constructemptymatrix(Kc, S)
  call fillmatrixfromtripletlist(Kc, tlist)
  call convertmatrixtoreal(Kc, K)

  if (isroot()) then
      write(*,*) "Patch Done"
  end if

  !! Print trace of patched matrix
  call dotmatrix(K, S, current_trace)
  if (isroot()) then
     write(*,*) "Current trace:", current_trace
     write(*,*) "Target trace:", target_trace
     write(*,*) "Trace error:", abs(current_trace - target_trace)
  end if

  !! Undo the permutation to the original ordering of the matrix
  call undopermutematrix(K, K, perm)

  !! Write to file
  ! call writematrixtomatrixmarket(K, o_file)

  !! Cleanup
  call destroy_negf(pnegf)
  call destructprocessgrid
  call mpi_finalize(ierr)
contains
  !> Wrap up setting up libNEGF
  subroutine setup_negf(pnegf)
    use mpi_f08
    type(Tnegf) :: pnegf
    call init_negf(pnegf)
    call set_mpi_bare_comm(pnegf, MPI_COMM_WORLD)
  end subroutine
  !> Read in the file with the ordering.
  subroutine get_order(idx_file, order, mid, lef1, lef2, rig1, rig2, mu)
    !> The name of the file containing the indexing information
    character(len=*), intent(in) :: idx_file
    !> The order of indices to permute the matrix to
    integer, dimension(:), allocatable, intent(out) :: order
    !> The last index in each region
    integer, intent(out) :: mid, lef1, lef2, rig1, rig2
    !> The chemical potential
    real(ntreal), intent(out) :: mu
    !! local variables
    integer, parameter :: fh = 16
    integer :: ii, length

    open(fh, file = idx_file, status = 'old')
    !! number of elements
    read(fh, *) length
    allocate(order(length))
    !! indices
    read(fh, *) order
    !! separators
    read(fh, *) mid, lef1, lef2, rig1, rig2
    !! chemical potential
    read(fh, *) mu
    close(fh)

  end subroutine
  !> Test if a row / column is in the appropriate region.
  logical function in_region(col, row)
    integer, intent(in) :: col, row
    !! LibNEGF computes: Device×Device + full Device×Contact_PL1 coupling
    in_region = (col <= mid .and. row <= mid) .or. &           ! Device×Device (1-48, 1-48)
        (col > mid .and. col <= lef1 .and. row <= mid) .or. &  ! Device×LEF:1 (full coupling)
        (row > mid .and. row <= lef1 .and. col <= mid) .or. &  ! LEF:1×Device (full coupling)
        (col > lef2 .and. col <= rig1 .and. row <= mid) .or. & ! Device×RIG:1 (full coupling)
        (row > lef2 .and. row <= rig1 .and. col <= mid)        ! RIG:1×Device (full coupling)
        ! (row <= rig1 .and. row > lef2 .and. col <= rig1 - lef2)
    ! in_region = (col <= mid  .and. row <= mid) .or. &
    !      (col <= lef1 .and. col > mid .and. row < lef1 - mid) .or. &
    !      (row <= lef1 .and. row > mid .and. col < lef1 - mid) .or. &
    !      (col <= rig1 .and. col > lef2 .and. row < rig1 - lef2) .or. &
    !      (row <= rig1 .and. row > lef2 .and. col < rig1 - lef2)
    ! in_region =  (col <= rig1 .and. col > lef2 .and. row < rig1 - lef2)
    ! in_region = &
    !      (col <= mid  .and. row <= mid)        .or. &
    !      (col <= lef1 .and. row >  lef1 - mid) .or. &
    !      (row <= lef1 .and. col >  lef1 - mid) .or. &
    !      (row > lef2  .and. row <= rig1 .and. &
    !       col > mid - (lef2 - mid) .and. col <= mid) .or. &
    !      (col > lef2  .and. col <= rig1 .and. &
    !       row > mid   - (lef2 - mid) .and. row <= mid)
  end function
  !> Extract all elements from libNEGF density (no filtering)
  subroutine extract_all_lnegf(ncol, rowpnt, colind, values)
    integer, intent(in)    :: ncol
    integer, intent(in)    :: rowpnt(:), colind(:)
    complex(kind=8), intent(in) :: values(:)
    integer :: jj, kk, row

    do jj = 1, ncol
      do kk = rowpnt(jj), rowpnt(jj+1) - 1
        row = colind(kk)
        trip%index_column = jj
        trip%index_row    = row
        trip%point_value  = values(kk)
        call appendtotripletlist(tlist, trip)
      end do
    end do
  end subroutine
  !> Extract elements from libNEGF density in specific regions
  subroutine extract_lnegf(ncol, rowpnt, colind, values)
    integer, intent(in)    :: ncol
    integer, intent(in)    :: rowpnt(:), colind(:)
    complex(kind=8), intent(in) :: values(:)
    integer :: jj, kk, row, count_extracted

    count_extracted = 0
    do jj = 1, ncol
      do kk = rowpnt(jj), rowpnt(jj+1) - 1
        row = colind(kk)
        if (in_region(jj, row) .eqv. .true.) then
          trip%index_column = jj
          trip%index_row    = row
          trip%point_value  = values(kk)
          call appendtotripletlist(tlist, trip)
          count_extracted = count_extracted + 1
        end if
      end do
    end do
    write(*,*) "  Extracted from libNEGF:", count_extracted, "elements"
  end subroutine
  subroutine extract_ntpoly(ncol, rowpnt, colind, values)
    integer, intent(in)    :: ncol
    integer, intent(in)    :: rowpnt(:), colind(:)
    complex(kind=8), intent(in) :: values(:)
    integer :: jj, kk, row, count_extracted

    count_extracted = 0
    do jj = 1, ncol
      !! NTPoly uses 0-based rowpnt, convert to 1-based for Fortran array access
      do kk = rowpnt(jj) + 1, rowpnt(jj+1)
        row = colind(kk)
        if (in_region(jj, row) .eqv. .false.) then
          trip%index_column = jj
          trip%index_row    = row
          trip%point_value  = values(kk)
          call appendtotripletlist(tlist, trip)
          count_extracted = count_extracted + 1
        end if
      end do
    end do
    write(*,*) "  Extracted from BigDFT:", count_extracted, "elements"
  end subroutine
  !> Compute trace of libNEGF density matrix
  subroutine check_lnegf_trace(ncol, rho_rowpnt, rho_colind, rho_values, &
                                s_ncol, s_rowpnt, s_colind, s_values)
    integer, intent(in) :: ncol, s_ncol
    integer, intent(in) :: rho_rowpnt(:), rho_colind(:)
    integer, intent(in) :: s_rowpnt(:), s_colind(:)
    complex(kind=8), intent(in) :: rho_values(:), s_values(:)
    complex(kind=8) :: trace
    integer :: jj, kk_rho, kk_s, row

    trace = (0.0d0, 0.0d0)
    !! Compute trace(rho * S) by summing diagonal elements
    do jj = 1, ncol
      do kk_rho = rho_rowpnt(jj), rho_rowpnt(jj+1) - 1
        row = rho_colind(kk_rho)
        !! Only count diagonal elements (weighted by S diagonal)
        if (row == jj) then
          !! Find S(jj, jj)
          do kk_s = s_rowpnt(jj) + 1, s_rowpnt(min(jj+1, s_ncol))
            if (s_colind(kk_s) == jj) then
              trace = trace + rho_values(kk_rho) * s_values(kk_s)
              exit
            end if
          end do
        end if
      end do
    end do
    write(*,*) "Trace of libNEGF density (diag only):", trace
  end subroutine
end program
