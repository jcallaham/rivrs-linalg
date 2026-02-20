! Driver to call SPRAL's SSIDS with user-supplied ordering and scaling,
! then dump factor statistics, D entries, and backward error for comparison
! with our Rust implementation.
!
! Input (stdin):
!   n nnz_lower
!   col_ptr[1] .. col_ptr[n+1]    (1-indexed, one per line)
!   row[1] val[1]                  (one per line, row 1-indexed)
!   ...
!   row[nnz_lower] val[nnz_lower]
!   ORDERING                       (marker line)
!   order[1] .. order[n]           (1-indexed positions, one per line)
!   SCALING                        (marker line)
!   scale[1] .. scale[n]           (scaling factors, one per line)
!   RHS                            (marker line)
!   rhs[1] .. rhs[n]               (right-hand side, one per line)
!
! Output (stdout):
!   INFORM section with factorization statistics
!   D section with D^{-1} entries (in pivot order)
!   SOLUTION section with solution vector
!   BACKWARD_ERROR with computed backward error
!
! Compile (after running tools/build_spral.sh):
!   METIS_LIB=$(find /workspace/rivrs-linalg/sparse/target -name "libmetis.a" | head -1)
!   gfortran -O2 -I /tmp/spral_ssids -o /tmp/spral_full_solve \
!     tools/spral_full_solve.f90 -L/tmp/spral_ssids -lspral \
!     $METIS_LIB -lopenblas -lstdc++ -lm

program spral_full_solve_driver
  use spral_ssids
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, k
  integer, allocatable :: ptr(:), row(:), order(:), piv_order(:)
  real(wp), allocatable :: val(:), scaling(:), rhs(:), x(:)
  real(wp), allocatable :: d(:,:)
  character(len=80) :: line

  type(ssids_akeep) :: akeep
  type(ssids_fkeep) :: fkeep
  type(ssids_options) :: options
  type(ssids_inform) :: inform

  ! Backward error computation
  real(wp), allocatable :: ax(:)
  real(wp) :: norm_a, norm_x, norm_b, norm_r, be
  integer :: j
  integer(kind=8) :: col_start, col_end

  ! Read dimensions
  read(*, *) n, nnz

  allocate(ptr(n+1), row(nnz), val(nnz))
  allocate(order(n), scaling(n), rhs(n), x(n))
  allocate(piv_order(n), d(2,n))
  allocate(ax(n))

  ! Read column pointers (1-indexed, lower triangle)
  do i = 1, n+1
    read(*, *) ptr(i)
  end do

  ! Read row indices and values
  do k = 1, nnz
    read(*, *) row(k), val(k)
  end do

  ! Read ordering marker
  read(*, '(a)') line
  if (trim(adjustl(line)) /= 'ORDERING') then
    write(0, '(a,a)') 'Expected ORDERING marker, got: ', trim(line)
    stop 1
  end if

  ! Read ordering (1-indexed: order(i) = position of variable i)
  do i = 1, n
    read(*, *) order(i)
  end do

  ! Read scaling marker
  read(*, '(a)') line
  if (trim(adjustl(line)) /= 'SCALING') then
    write(0, '(a,a)') 'Expected SCALING marker, got: ', trim(line)
    stop 1
  end if

  ! Read scaling factors (original variable order)
  do i = 1, n
    read(*, *) scaling(i)
  end do

  ! Read RHS marker
  read(*, '(a)') line
  if (trim(adjustl(line)) /= 'RHS') then
    write(0, '(a,a)') 'Expected RHS marker, got: ', trim(line)
    stop 1
  end if

  ! Read right-hand side
  do i = 1, n
    read(*, *) rhs(i)
  end do

  write(0, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Configure SPRAL SSIDS
  options%ordering = 0       ! user-supplied ordering
  options%scaling = 0        ! user-supplied scaling (via scale argument)
  options%print_level = 0    ! minimal output
  options%action = .true.    ! continue on singular

  ! Analyse phase
  call ssids_analyse(.true., n, ptr, row, akeep, options, inform, order=order)
  if (inform%flag < 0) then
    write(0, '(a,i5)') 'ssids_analyse error flag=', inform%flag
    write(*, '(a,i5)') 'ANALYSE_FLAG ', inform%flag
    stop 1
  end if
  write(0, '(a,i5)') 'analyse flag=', inform%flag

  ! Factor phase (indefinite)
  call ssids_factor(.false., val, akeep, fkeep, options, inform, scale=scaling)
  if (inform%flag < 0) then
    write(0, '(a,i5)') 'ssids_factor error flag=', inform%flag
    write(*, '(a,i5)') 'FACTOR_FLAG ', inform%flag
    stop 1
  end if
  write(0, '(a,i5)') 'factor flag=', inform%flag

  ! Dump inform statistics
  write(*, '(a)') 'INFORM'
  write(*, '(a,i8)') 'flag ', inform%flag
  write(*, '(a,i8)') 'matrix_rank ', inform%matrix_rank
  write(*, '(a,i8)') 'num_delay ', inform%num_delay
  write(*, '(a,i8)') 'num_two ', inform%num_two
  write(*, '(a,i8)') 'num_neg ', inform%num_neg
  write(*, '(a,i8)') 'maxfront ', inform%maxfront
  write(*, '(a,i8)') 'maxsupernode ', inform%maxsupernode
  write(*, '(a,i8)') 'num_sup ', inform%num_sup
  write(*, '(a,i20)') 'num_factor ', inform%num_factor
  write(*, '(a,i20)') 'num_flops ', inform%num_flops
  write(*, '(a,i8)') 'not_first_pass ', inform%not_first_pass
  write(*, '(a,i8)') 'not_second_pass ', inform%not_second_pass

  ! Enquire D^{-1} and pivot order
  call ssids_enquire_indef(akeep, fkeep, options, inform, &
       piv_order=piv_order, d=d)

  write(*, '(a)') 'PIVOT_ORDER'
  do i = 1, n
    write(*, '(i8)') piv_order(i)
  end do

  write(*, '(a)') 'D_INV'
  do i = 1, n
    write(*, '(es24.16,1x,es24.16)') d(1,i), d(2,i)
  end do

  ! Solve phase
  x(1:n) = rhs(1:n)
  call ssids_solve(x, akeep, fkeep, options, inform)
  if (inform%flag < 0) then
    write(0, '(a,i5)') 'ssids_solve error flag=', inform%flag
    write(*, '(a,i5)') 'SOLVE_FLAG ', inform%flag
    stop 1
  end if

  write(*, '(a)') 'SOLUTION'
  do i = 1, n
    write(*, '(es24.16)') x(i)
  end do

  ! Compute backward error: ||Ax - b|| / (||A|| ||x|| + ||b||)
  ! A is stored as lower triangle CSC, but represents full symmetric matrix
  ax(1:n) = 0.0_wp
  norm_a = 0.0_wp
  do j = 1, n
    col_start = ptr(j)
    col_end = ptr(j+1) - 1
    do k = int(col_start), int(col_end)
      i = row(k)
      ! Lower triangle: a(i,j)
      ax(i) = ax(i) + val(k) * x(j)
      norm_a = max(norm_a, abs(val(k)))
      if (i /= j) then
        ! Upper triangle (symmetric): a(j,i)
        ax(j) = ax(j) + val(k) * x(i)
      end if
    end do
  end do

  ! Residual r = Ax - b
  norm_r = 0.0_wp
  norm_x = 0.0_wp
  norm_b = 0.0_wp
  do i = 1, n
    norm_r = max(norm_r, abs(ax(i) - rhs(i)))
    norm_x = max(norm_x, abs(x(i)))
    norm_b = max(norm_b, abs(rhs(i)))
  end do

  if (norm_a * norm_x + norm_b > 0.0_wp) then
    be = norm_r / (norm_a * norm_x + norm_b)
  else
    be = 0.0_wp
  end if

  write(*, '(a,es24.16)') 'BACKWARD_ERROR ', be

  write(0, '(a,es12.4)') 'Backward error: ', be
  write(0, '(a,i8,a,i8,a,i8,a,i8)') &
    'delays=', inform%num_delay, &
    ' 2x2=', inform%num_two, &
    ' not_1st=', inform%not_first_pass, &
    ' not_2nd=', inform%not_second_pass

  ! Cleanup
  call ssids_free(akeep, fkeep, i)

end program spral_full_solve_driver
