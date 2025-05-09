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
  real(ntreal) :: mu
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

  !! Setup MPI, NTPoly, libNEGF
  call mpi_init_thread(mpi_thread_serialized, prov, ierr)
  call constructprocessgrid(mpi_comm_world)
  call setup_negf(pnegf)

  !! Read in the command line arguments
  call get_command_argument(1, h_file)
  call get_command_argument(2, s_file)
  call get_command_argument(3, k_file)
  call get_command_argument(4, idx_file)
  call get_command_argument(5, temps)
  read(temps, *) mu
  call get_command_argument(6, o_file)

  !! Read in the BigDFT Matrices
  call constructmatrixfrommatrixmarket(H, h_file)
  call constructmatrixfrommatrixmarket(S, s_file)
  call constructmatrixfrommatrixmarket(K, k_file)

  !! The number of occupied orbitals read from the trace
  call dotmatrix(K, S, dotv)
  if (isroot()) then
      write(*,*) "Trace should be", dotv
  end if

  !! Read in the indexing information
  call get_order(idx_file, order, mid, lef1, lef2, rig1, rig2)
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
     perm%reverse_index_lookup(perm%index_lookup(II)) = II 
  end do
  call permutematrix(H, H, perm)
  call permutematrix(S, S, perm)
  call permutematrix(K, K, perm)

  !! snap to sparsity pattern
  call snapmatrixtosparsitypattern(H, K)
  call snapmatrixtosparsitypattern(S, K)

  !! convert to complex
  call convertmatrixtocomplex(H, Hc)
  call convertmatrixtocomplex(S, Sc)
  call convertmatrixtocomplex(K, Kc)

  !! Gather locally
  call gathermatrixtoprocess(Hc, Hloc)
  call gathermatrixtoprocess(Sc, Sloc)
  call gathermatrixtoprocess(Kc, Kloc)

  if (isroot()) then
      write(*,*) "Setup Done"
  end if

  !! Convert to the CSR convention of libNEGF
  Hloc%outer_index = Hloc%outer_index + 1
  Sloc%outer_index = Sloc%outer_index + 1

  !! Convert to libNEGF matrices
  call create_HS(pnegf, 1) ! 1 k-point
  call create_DM(pnegf, 1)
  call set_h(pnegf, Hloc%rows, Hloc%values, Hloc%inner_index, Hloc%outer_index)
  call set_s(pnegf, Sloc%rows, Sloc%values, Sloc%inner_index, Sloc%outer_index)

  !! Call the libNEGF solver
  call init_contacts(pnegf, 2)
  call init_structure(pnegf, 2, surfstart, surfend, contend, 1, plend, cblk)
  call get_params(pnegf, params)
  params%mu(1:2) = mu
  params%Emin = mu - 0.2
  params%Emax = mu + 0.2
  params%Estep = 2.d-2 
  params%kbT_dm(1:2) = 1e-3
  params%kbT_t(1:2) = 1e-3
  params%ec = -2.0
  call set_params(pnegf, params)
  call compute_density_dft(pnegf)
  call mpi_allreduce(mpi_in_place, pnegf%rho%nzval, &
       & size(pnegf%rho%nzval), MPI_DOUBLE_COMPLEX,&
       & MPI_SUM, MPI_COMM_WORLD, ierr)

  if (isroot()) then
      write(*,*) "Density Done"
  end if

  !! Compute the transmission
  call compute_current(pnegf) 
  call associate_transmission(pnegf, transmission) 
  call mpi_allreduce(mpi_in_place, transmission, &
       & size(transmission, 1) * size(transmission, 2), &
       & MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  if (isroot()) then
    call write_tunneling_and_dos(pnegf)
  end if

  if (isroot()) then
      write(*,*) "Transmission Done"
  end if

  !! Patch back in the libNEGF density update
  call constructtripletlist(tlist)
  if (isroot()) then
    call extract(pnegf%rho%ncol, pnegf%rho%rowpnt,&
                 pnegf%rho%colind, pnegf%rho%nzval, .true.)

    call extract(kloc%columns, kloc%outer_index, &
                 kloc%inner_index, kloc%values, .false.)
  end if
  call constructemptymatrix(Kc, S)
  call fillmatrixfromtripletlist(Kc, tlist)
  call convertmatrixtoreal(Kc, K)

  if (isroot()) then
      write(*,*) "Patch Done"
  end if

  !! Print trace of patched matrix to check if it is reasonable
  call dotmatrix(K, S, dotv)
  if (isroot()) then
     write(*,*) "Trace is:", dotv
  end if

  !! Undo the permutation to the original ordering of the matrix
  call undopermutematrix(K, K, perm)

  !! Write to file
  call writematrixtomatrixmarket(K, o_file)

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
  subroutine get_order(idx_file, order, mid, lef1, lef2, rig1, rig2)
    !> The name of the file containing the indexing information
    character(len=*), intent(in) :: idx_file
    !> The order of indices to permute the matrix to
    integer, dimension(:), allocatable, intent(out) :: order
    !> The last index in each region
    integer, intent(out) :: mid, lef1, lef2, rig1, rig2
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
    close(fh)

  end subroutine
  !> Test if a row / column is in the appropriate region.
  logical function in_region(col,row)
    integer, intent(in) :: col, row
    in_region = &
         (col<=mid .and. row<=mid)        .or. &
         (col<=lef1 .and. row> lef1-mid) .or. &
         (row<=lef1 .and. col> lef1-mid) .or. &
         (row> lef2 .and. row<=rig1 .and. &
          col> mid-(lef2-mid) .and. col<=mid) .or. &
         (col> lef2 .and. col<=rig1 .and. &
          row> mid-(lef2-mid) .and. row<=mid)
  end function
  !> 
  subroutine extract(ncol, rowpnt, colind, values, want_inside)
    integer, intent(in)    :: ncol
    integer, intent(in)    :: rowpnt(:), colind(:)
    complex(kind=8), intent(in) :: values(:)
    logical, intent(in)    :: want_inside
    integer :: jj, kk, row

    do jj = 1, ncol
      do kk = rowpnt(jj)+1, rowpnt(jj+1) - 1
        row = colind(kk)
        if (in_region(jj, row) .eqv. want_inside) then
          trip%index_column = jj
          trip%index_row    = row
          trip%point_value  = values(kk)
          call appendtotripletlist(tlist, trip)
        end if
      end do
    end do
  end subroutine
end program
