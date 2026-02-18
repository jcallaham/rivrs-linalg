! Driver to dump SPRAL's scaling factors and duals for comparison with our Rust code.
!
! Compile:
!   gfortran -O2 -o spral_scaling_dump \
!     /opt/references/spral/src/matrix_util.f90 \
!     /opt/references/spral/src/scaling.f90 \
!     tools/spral_scaling_dump.f90

program spral_scaling_dump
  use spral_scaling
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, j, k
  integer, allocatable :: ptr(:), row(:), match_perm(:)
  real(wp), allocatable :: val(:), scaling(:)
  real(wp) :: v, max_violation, max_scaled

  type(hungarian_options) :: options
  type(hungarian_inform) :: inform

  ! Read dimensions
  read(*, *) n, nnz

  allocate(ptr(n+1), row(nnz), val(nnz))
  allocate(scaling(n), match_perm(n))

  ! Read column pointers (1-indexed)
  do i = 1, n+1
    read(*, *) ptr(i)
  end do

  ! Read row indices and values
  do k = 1, nnz
    read(*, *) row(k), val(k)
  end do

  write(0, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Call SPRAL Hungarian scaling
  call hungarian_scale_sym(n, ptr, row, val, scaling, options, inform, match=match_perm)

  write(0, '(a,i5)') 'inform%flag = ', inform%flag
  write(0, '(a,i8)') 'inform%matched = ', inform%matched

  if (inform%flag < 0) then
    write(0, '(a)') 'SPRAL returned error!'
    stop 1
  end if

  ! Dump scaling factors (one per line)
  write(*, '(a,i8)') 'N ', n
  write(*, '(a,i8)') 'MATCHED ', inform%matched
  do i = 1, n
    write(*, '(es24.16)') scaling(i)
  end do

  ! Dump matching (1-indexed)
  write(*, '(a)') 'MATCHING'
  do i = 1, n
    write(*, '(i8)') match_perm(i)
  end do

  ! Check violations
  max_violation = 0.0d0
  do j = 1, n
    do k = ptr(j), ptr(j+1)-1
      i = row(k)
      v = abs(scaling(i) * val(k) * scaling(j))
      if (v - 1.0d0 > max_violation) max_violation = v - 1.0d0
    end do
  end do
  write(0, '(a,es12.4)') 'max_violation = ', max_violation

end program spral_scaling_dump
