! Minimal driver to test SPRAL's hungarian_scale_sym on matrices
! exported from our Rust code.
!
! Input format (text file on stdin):
!   Line 1: n nnz_lower
!   Next n+1 lines: ptr(1:n+1)  (1-indexed column pointers for lower triangle)
!   Next nnz_lower lines: row(1:nnz_lower) val(1:nnz_lower)  (row index + value)
!
! Compile:
!   gfortran -O2 -o spral_scaling_test \
!     /opt/references/spral/src/matrix_util.f90 \
!     /opt/references/spral/src/scaling.f90 \
!     tools/spral_scaling_test.f90

program spral_scaling_test
  use spral_scaling
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, j, k, info_stat
  integer, allocatable :: ptr(:), row(:)
  real(wp), allocatable :: val(:), scaling(:), rmax(:)
  integer, allocatable :: match(:)
  real(wp) :: v, max_violation, err_tol
  integer :: violation_count

  type(hungarian_options) :: options
  type(hungarian_inform) :: inform

  err_tol = 5.0d-14

  ! Read dimensions
  read(*, *) n, nnz

  allocate(ptr(n+1), row(nnz), val(nnz))
  allocate(scaling(n), match(n), rmax(n))

  ! Read column pointers (1-indexed)
  do i = 1, n+1
    read(*, *) ptr(i)
  end do

  ! Read row indices and values
  do k = 1, nnz
    read(*, *) row(k), val(k)
  end do

  write(*, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Call SPRAL Hungarian scaling
  call hungarian_scale_sym(n, ptr, row, val, scaling, options, inform, match=match)

  write(*, '(a,i5)') 'inform%flag = ', inform%flag
  if (inform%flag < 0) then
    write(*, '(a)') 'SPRAL returned error!'
    stop 1
  end if
  write(*, '(a,i8)') 'inform%matched = ', inform%matched

  ! Check scaling bound: |s_i * a_ij * s_j| <= 1 + tol
  max_violation = 0.0d0
  violation_count = 0
  rmax = 0.0d0

  do j = 1, n
    do k = ptr(j), ptr(j+1)-1
      i = row(k)
      v = abs(scaling(i) * val(k) * scaling(j))
      ! Track row max (lower triangle contributes to both i and j)
      if (v > rmax(i)) rmax(i) = v
      if (i /= j .and. v > rmax(j)) rmax(j) = v
      if (v > 1.0d0 + err_tol) then
        violation_count = violation_count + 1
        if (v - 1.0d0 > max_violation) max_violation = v - 1.0d0
        if (violation_count <= 5) then
          write(*, '(a,i8,a,i8,a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
            '  violation: (', i, ',', j, ') scaled=', v, &
            ' s[i]=', scaling(i), ' s[j]=', scaling(j), ' val=', val(k)
        end if
      end if
    end do
  end do

  write(*, '(a,i8)') 'Total violations (|s*a*s| > 1+5e-14): ', violation_count
  write(*, '(a,es12.4)') 'Max violation excess: ', max_violation

  ! Check row max quality
  do i = 1, n
    if (rmax(i) > 0.0d0 .and. rmax(i) < 1.0d0 - err_tol) then
      write(*, '(a,i8,a,es12.4)') '  rmax violation: row ', i, ' rmax=', rmax(i)
    end if
  end do

  write(*, '(a)') 'Done.'

end program spral_scaling_test
