! Dump SPRAL's raw dual variables (dualu, dualv) after Hungarian matching.
!
! This replicates hungarian_wrapper's preprocessing and calls hungarian_match
! directly so we can access the raw duals before they get combined into scaling.
!
! Compile:
!   gfortran -O2 -o spral_dump_duals \
!     /opt/references/spral/src/matrix_util.f90 \
!     /opt/references/spral/src/scaling.f90 \
!     tools/spral_dump_duals.f90

program spral_dump_duals
  use spral_matrix_util, only : half_to_full
  use spral_scaling, only : hungarian_match
  implicit none

  integer, parameter :: wp = kind(0d0)
  real(wp), parameter :: zero = 0.0_wp
  integer :: n, nnz_lower, i, j, matched
  integer(kind=8) :: klong, jlong, ne
  integer, allocatable :: ptr_lower(:), row_lower(:)
  real(wp), allocatable :: val_lower(:)
  integer(kind=8), allocatable :: ptr2(:)
  integer, allocatable :: row2(:), iw(:), match_perm(:)
  real(wp), allocatable :: val2(:), dualu(:), dualv(:), cmax(:)
  real(wp) :: colmax, v
  integer :: stat

  ! Read dimensions
  read(*, *) n, nnz_lower

  allocate(ptr_lower(n+1), row_lower(nnz_lower), val_lower(nnz_lower))

  ! Read column pointers (1-indexed)
  do i = 1, n+1
    read(*, *) ptr_lower(i)
  end do

  ! Read row indices and values
  do i = 1, nnz_lower
    read(*, *) row_lower(i), val_lower(i)
  end do

  write(0, '(a,i8,a,i10)') 'Lower triangle: n=', n, ' nnz=', nnz_lower

  ! Replicate hungarian_wrapper preprocessing (lines 623-657)
  ne = 2 * (ptr_lower(n+1) - 1)
  allocate(ptr2(n+1), row2(ne), val2(ne), iw(5*n))
  allocate(dualu(n), dualv(n), cmax(n), match_perm(n))

  ! Copy and take log of abs values (lines 636-648)
  klong = 1
  do i = 1, n
    ptr2(i) = klong
    do jlong = ptr_lower(i), ptr_lower(i+1)-1
      if (val_lower(jlong) .eq. zero) cycle
      row2(klong) = row_lower(jlong)
      val2(klong) = abs(val_lower(jlong))
      klong = klong + 1
    end do
    val2(ptr2(i):klong-1) = log(val2(ptr2(i):klong-1))
  end do
  ptr2(n+1) = klong

  ! Expand to full symmetric (line 650)
  call half_to_full(n, row2, ptr2, iw, a=val2)

  write(0, '(a,i10)') 'Full nnz after expansion: ', ptr2(n+1)-1

  ! Compute column maxima and costs (lines 653-657)
  do i = 1, n
    colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
    cmax(i) = colmax
    val2(ptr2(i):ptr2(i+1)-1) = colmax - val2(ptr2(i):ptr2(i+1)-1)
  end do

  ! Call hungarian_match directly (line 659)
  call hungarian_match(n, n, ptr2, row2, val2, match_perm, matched, dualu, &
       dualv, stat)

  write(0, '(a,i8)') 'matched = ', matched
  write(0, '(a,i5)') 'stat = ', stat

  ! Check dual feasibility: u[i] + v[j] <= c[i,j] + eps for all (i,j)
  ! where c[i,j] is stored in val2
  v = 0.0_wp
  do j = 1, n
    do klong = ptr2(j), ptr2(j+1)-1
      i = row2(klong)
      ! Reduced cost should be >= 0 for feasibility
      ! c[i,j] - u[i] - v[j] >= 0
      ! i.e., u[i] + v[j] <= c[i,j]
      v = max(v, dualu(i) + dualv(j) - val2(klong))
    end do
  end do
  write(0, '(a,es12.4)') 'max dual infeasibility (u+v-c): ', v

  ! Dump duals
  write(*, '(a,i8)') 'N ', n
  write(*, '(a,i8)') 'MATCHED ', matched

  write(*, '(a)') 'DUALU'
  do i = 1, n
    write(*, '(es24.16)') dualu(i)
  end do

  write(*, '(a)') 'DUALV'
  do i = 1, n
    write(*, '(es24.16)') dualv(i)
  end do

  write(*, '(a)') 'CMAX'
  do i = 1, n
    write(*, '(es24.16)') cmax(i)
  end do

  write(*, '(a)') 'MATCHING'
  do i = 1, n
    write(*, '(i8)') match_perm(i)
  end do

end program spral_dump_duals
