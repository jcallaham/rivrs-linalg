! Driver to call SPRAL's match_order_metis() and dump ordering + scaling
! for comparison with our Rust implementation.
!
! Input (stdin): lower-triangle CSC format
!   n nnz_lower
!   col_ptr[1] .. col_ptr[n+1]    (1-indexed, one per line)
!   row[1] val[1]                  (one per line, row 1-indexed)
!   ...
!   row[nnz_lower] val[nnz_lower]
!
! Output (stdout):
!   N <n>
!   FLAG <flag>
!   SCALING                     (followed by n lines of scaling factors)
!   ORDERING                    (followed by n lines of elimination order)
!
! Compile:
!   # First compile SPRAL modules into /tmp/spral_mo/:
!   gfortran -O2 -c -o /tmp/spral_mo/matrix_util.o /opt/references/spral/src/matrix_util.f90 -J /tmp/spral_mo
!   gfortran -O2 -c -o /tmp/spral_mo/scaling.o /opt/references/spral/src/scaling.f90 -J /tmp/spral_mo -I /tmp/spral_mo
!   gfortran -O2 -c -o /tmp/spral_mo/metis5_wrapper.o -DSPRAL_HAVE_METIS_H=0 /opt/references/spral/src/metis5_wrapper.F90 -J /tmp/spral_mo -I /tmp/spral_mo
!   gfortran -O2 -c -o /tmp/spral_mo/match_order.o /opt/references/spral/src/match_order.f90 -J /tmp/spral_mo -I /tmp/spral_mo
!   # Then compile and link:
!   METIS_LIB=$(find /workspace/rivrs-linalg/sparse/target -name "libmetis.a" | head -1)
!   gfortran -O2 -I /tmp/spral_mo -o /tmp/spral_match_order \
!     /tmp/spral_mo/matrix_util.o /tmp/spral_mo/scaling.o /tmp/spral_mo/metis5_wrapper.o /tmp/spral_mo/match_order.o \
!     tools/spral_match_order.f90 $METIS_LIB -lm

program spral_match_order_driver
  use spral_match_order
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, k
  integer, allocatable :: ptr(:), row(:), order(:)
  real(wp), allocatable :: val(:), scaling(:)
  integer :: flag, stat

  ! Read dimensions
  read(*, *) n, nnz

  allocate(ptr(n+1), row(nnz), val(nnz))
  allocate(scaling(n), order(n))

  ! Read column pointers (1-indexed)
  do i = 1, n+1
    read(*, *) ptr(i)
  end do

  ! Read row indices and values
  do k = 1, nnz
    read(*, *) row(k), val(k)
  end do

  write(0, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Call SPRAL match_order_metis (expects full matrix, both triangles)
  call match_order_metis(n, ptr, row, val, order, scaling, flag, stat)

  write(0, '(a,i5)') 'flag = ', flag
  write(0, '(a,i5)') 'stat = ', stat

  ! Dump results
  write(*, '(a,i8)') 'N ', n
  write(*, '(a,i5)') 'FLAG ', flag

  if (flag < 0) then
    write(0, '(a)') 'SPRAL match_order_metis returned error!'
    stop 1
  end if

  ! Scaling factors (linear domain — match_order_metis already exp()'s them)
  write(*, '(a)') 'SCALING'
  do i = 1, n
    write(*, '(es24.16)') scaling(i)
  end do

  ! Elimination ordering: order(i) = position of variable i
  write(*, '(a)') 'ORDERING'
  do i = 1, n
    write(*, '(i8)') order(i)
  end do

end program spral_match_order_driver
