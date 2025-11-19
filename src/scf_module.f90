!> Module for self-consistent chemical potential search using NEGF
!> Provides clean interface for SCF calculations with NTPoly matrices
module scf_module
  !! libNEGF includes
  use libnegf
  use lib_param
  !! NTPoly includes
  use datatypesmodule
  use matrixconversionmodule
  use permutationmodule
  use processgridmodule
  use psmatrixalgebramodule
  use psmatrixmodule
  use smatrixalgebramodule
  use smatrixmodule
  use tripletlistmodule
  use tripletmodule
  !! MPI
  use mpi
  implicit none
  private

  !! Public types
  public :: scf_params_t, region_info_t
  !! Public routines
  public :: scf_solve

  !> SCF convergence parameters
  type :: scf_params_t
    integer :: max_iterations = 15
    real(ntreal) :: tolerance = 1.0d-3
    real(ntreal) :: bracket_width = 0.05d0
    real(ntreal) :: mu_initial = 0.0d0
  end type scf_params_t

  !> Region boundary information for libNEGF
  type :: region_info_t
    integer :: mid    ! Last index of central region
    integer :: lef1   ! Last index of left buffer
    integer :: lef2   ! Last index of left contact
    integer :: rig1   ! Last index of right buffer
    integer :: rig2   ! Last index of right contact
  end type region_info_t

contains

  !> Main SCF solver routine
  !> Takes distributed NTPoly matrices and returns converged density
  subroutine scf_solve(H, S, K, regions, params, mu_final, converged, num_iter)
    type(matrix_ps), intent(in) :: H, S
    type(matrix_ps), intent(inout) :: K
    type(region_info_t), intent(in) :: regions
    type(scf_params_t), intent(in) :: params
    real(ntreal), intent(out) :: mu_final
    logical, intent(out) :: converged
    integer, intent(out) :: num_iter

    !! Local SCF variables
    real(ntreal) :: mu, mu_low, mu_high
    real(ntreal) :: trace_low, trace_high
    real(ntreal) :: target_trace, current_trace, trace_error
    logical :: have_trace_low, have_trace_high
    integer :: iteration

    !! Pre-converted matrices (done ONCE outside loop)
    type(matrix_ps) :: Hc, Sc, Kc
    type(matrix_lsc) :: Hloc, Sloc, Kloc_contacts

    !! Pre-built matrix with contact regions (done ONCE, never change!)
    type(matrix_ps) :: K_contacts
    type(tripletlist_c) :: contact_tlist
    type(triplet_c) :: trip

    !! Initialize bracket around DFT mu
    mu_low = params%mu_initial - params%bracket_width
    mu_high = params%mu_initial + params%bracket_width
    have_trace_low = .false.
    have_trace_high = .false.

    !! Compute target trace from original K and S
    call dotmatrix(K, S, target_trace)

    !! ========================================================================
    !! ONE-TIME SETUP: Convert and gather H, S (they never change!)
    !! ========================================================================
    if (isroot()) write(*,*) "Converting H and S to complex and gathering (ONCE)..."
    call convertmatrixtocomplex(H, Hc)
    call convertmatrixtocomplex(S, Sc)
    call gathermatrixtoprocess(Hc, Hloc)
    call gathermatrixtoprocess(Sc, Sloc)

    !! ========================================================================
    !! ONE-TIME SETUP: Build K_contacts matrix (lef2, rig2 - never change!)
    !! ========================================================================
    if (isroot()) write(*,*) "Building K_contacts matrix with unchanging regions (ONCE)..."
    call convertmatrixtocomplex(K, Kc)
    call gathermatrixtoprocess(Kc, Kloc_contacts)

    !! Extract contact regions and build K_contacts matrix
    call constructtripletlist(contact_tlist)
    if (isroot()) then
      call extract_ntpoly(Kloc_contacts%columns, Kloc_contacts%outer_index, &
                          Kloc_contacts%inner_index, Kloc_contacts%values, &
                          regions%mid, regions%lef1, regions%lef2, &
                          regions%rig1, regions%rig2, contact_tlist, trip)
    end if

    !! Build distributed matrix from contact regions
    call convertmatrixtoreal(Kc, K_contacts)  ! Start with template
    call constructemptymatrix(Kc, S)
    call fillmatrixfromtripletlist(Kc, contact_tlist)
    call convertmatrixtoreal(Kc, K_contacts)

    if (isroot()) then
      write(*,*)
      write(*,*) "=== Starting SCF Iterations ==="
      write(*,'(A,F12.6)') " Target trace: ", target_trace
      write(*,'(A,F12.8,A)') " Initial mu: ", params%mu_initial, " Ha"
      write(*,'(A,F10.6,A,F10.6,A)') " Bracket: [", mu_low, ", ", mu_high, "]"
      write(*,*)
    end if

    !! SCF loop
    converged = .false.
    do iteration = 1, params%max_iterations
      num_iter = iteration

      !! Choose next mu using linear interpolation or bisection
      if (have_trace_low .and. have_trace_high) then
        !! Linear interpolation (regula falsi)
        if (abs(trace_high - trace_low) > 1.0d-12) then
          mu = mu_low + (mu_high - mu_low) * (target_trace - trace_low) / (trace_high - trace_low)
          mu = max(mu_low, min(mu_high, mu))
        else
          mu = (mu_low + mu_high) / 2.0d0
        end if
      else
        !! Bisection until we have both bracket traces
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

      !! Run NEGF iteration (H, S pre-converted; K_contacts pre-built)
      call run_negf_iteration(Hloc, Sloc, K_contacts, S, K, regions, mu, current_trace)

      !! Check convergence
      trace_error = abs(current_trace - target_trace)

      if (isroot()) then
        write(*,'(A,F12.6,A,F12.6,A,F12.6,A)') " Trace: ", current_trace, &
             " (target: ", target_trace, ", error: ", trace_error, ")"
      end if

      if (trace_error < params%tolerance) then
        converged = .true.
        mu_final = mu
        if (isroot()) then
          write(*,*)
          write(*,*) "*** CONVERGED in", iteration, "iterations ***"
          write(*,'(A,F12.8,A)') " Final mu = ", mu, " Ha"
        end if
        exit
      end if

      !! Update bracket
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

    if (.not. converged) then
      mu_final = mu
      if (isroot()) then
        write(*,*)
        write(*,*) "*** Max iterations reached ***"
        write(*,'(A,F12.8,A,F12.8)') " Final mu = ", mu, " Ha, error = ", trace_error
      end if
    end if

  end subroutine scf_solve

  !> Run a single NEGF iteration and return updated density
  !> Hloc, Sloc are pre-converted to local complex format (done ONCE)
  !> K_contacts is a pre-built matrix with unchanging contact regions (done ONCE)
  !> We build K_device each iteration, then add: K = K_contacts + K_device
  subroutine run_negf_iteration(Hloc, Sloc, K_contacts, S, K, regions, mu_val, current_trace)
    type(matrix_lsc), intent(in) :: Hloc, Sloc
    type(matrix_ps), intent(in) :: K_contacts, S
    type(matrix_ps), intent(inout) :: K
    type(region_info_t), intent(in) :: regions
    real(ntreal), intent(in) :: mu_val
    real(ntreal), intent(out) :: current_trace

    !! Local variables - matrices
    type(matrix_ps) :: Kc, K_device

    !! libNEGF variables
    type(Tnegf), target :: pnegf
    type(lnParams) :: params
    real(kind(1.d0)), dimension(:,:), pointer :: transmission

    !! Region setup for libNEGF
    integer, dimension(2) :: surfstart, surfend, contend
    integer, dimension(1) :: plend
    integer, dimension(:), allocatable :: cblk

    !! Temporary variables
    integer :: ierr
    type(tripletlist_c) :: tlist
    type(triplet_c) :: trip

    !! Setup libNEGF
    call setup_negf_local(pnegf)

    !! Setup libNEGF structure
    allocate(cblk(2))
    plend(1) = regions%mid
    surfstart(1) = regions%mid + 1
    surfstart(2) = regions%lef2 + 1
    surfend(1) = regions%mid
    surfend(2) = regions%lef2
    contend(1) = regions%lef2
    contend(2) = regions%rig2
    cblk = 1

    !! Set H and S matrices (pre-converted, just setting pointers)
    call create_HS(pnegf, 1)
    call create_DM(pnegf, 1)
    call set_h(pnegf, Hloc%rows, Hloc%values, Hloc%inner_index, Hloc%outer_index + 1)
    call set_s(pnegf, Sloc%rows, Sloc%values, Sloc%inner_index, Sloc%outer_index + 1)

    !! Initialize contacts and structure
    call init_contacts(pnegf, 2)
    call init_structure(pnegf, 2, surfstart, surfend, contend, 1, plend, cblk)

    !! Set parameters with current mu
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
         size(pnegf%rho%nzval), MPI_DOUBLE_COMPLEX, &
         MPI_SUM, MPI_COMM_WORLD, ierr)

    !! Compute transmission
    call compute_current(pnegf)
    call associate_transmission(pnegf, transmission)
    call mpi_allreduce(mpi_in_place, transmission, &
         size(transmission, 1) * size(transmission, 2), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    !! Build K_device matrix with only libNEGF device regions (mid, lef1, rig1)
    call constructtripletlist(tlist)
    if (isroot()) then
      call extract_lnegf(pnegf%rho%ncol, pnegf%rho%rowpnt, &
                         pnegf%rho%colind, pnegf%rho%nzval, &
                         regions%mid, regions%lef1, regions%lef2, &
                         regions%rig1, regions%rig2, tlist, trip)
    end if

    call constructemptymatrix(Kc, S)
    call fillmatrixfromtripletlist(Kc, tlist)
    call convertmatrixtoreal(Kc, K_device)

    !! Add K_contacts + K_device to get full density matrix
    call incrementmatrix(K_contacts, K_device)
    call copymatrix(K_device, K)

    !! Compute current trace
    call dotmatrix(K, S, current_trace)

    !! Cleanup
    deallocate(cblk)
    call destroy_negf(pnegf)

  end subroutine run_negf_iteration

  !> Extract elements from libNEGF density in device regions
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

  !> Extract elements from NTPoly density in contact regions
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

  !> Test if matrix element is in libNEGF-computed region
  logical function in_region(col, row, mid, lef1, lef2, rig1, rig2)
    integer, intent(in) :: col, row, mid, lef1, lef2, rig1, rig2
    in_region = (col <= mid .and. row <= mid) .or. &
        (col > mid .and. col <= lef1 .and. row <= mid) .or. &
        (row > mid .and. row <= lef1 .and. col <= mid) .or. &
        (col > lef2 .and. col <= rig1 .and. row <= mid) .or. &
        (row > lef2 .and. row <= rig1 .and. col <= mid)
  end function in_region

  !> Setup libNEGF with MPI interface
  subroutine setup_negf_local(pnegf)
    use mpi_f08
    type(Tnegf) :: pnegf
    call init_negf(pnegf)
    call set_mpi_bare_comm(pnegf, MPI_COMM_WORLD)
  end subroutine setup_negf_local

end module scf_module
