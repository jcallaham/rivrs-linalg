! Run SPRAL's hungarian_match on a pre-built cost graph.
! This tests whether our cost graph construction is the issue.
!
! Input format (same as export_cost_graph.rs output):
!   Line 1: n nnz
!   Next n+1 lines: 1-indexed column pointers (int64)
!   Next nnz lines: 1-indexed row_index cost
!
! Compile:
!   gfortran -O2 -o spral_match_costgraph \
!     /opt/references/spral/src/matrix_util.f90 \
!     /opt/references/spral/src/scaling.f90 \
!     tools/spral_match_costgraph.f90

program spral_match_costgraph
  use spral_scaling, only : hungarian_match
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, j, matched, stat
  integer(kind=8) :: k
  integer(kind=8), allocatable :: ptr(:)
  integer, allocatable :: row_idx(:), match_perm(:)
  real(wp), allocatable :: cost(:), dualu(:), dualv(:)
  real(wp) :: v, max_infeasibility

  ! Read dimensions
  read(*, *) n, nnz

  allocate(ptr(n+1), row_idx(nnz), cost(nnz))
  allocate(dualu(n), dualv(n), match_perm(n))

  ! Read column pointers (1-indexed, as int64)
  do i = 1, n+1
    read(*, *) ptr(i)
  end do

  ! Read row indices and costs
  do k = 1, nnz
    read(*, *) row_idx(k), cost(k)
  end do

  write(0, '(a,i8,a,i10)') 'Cost graph: n=', n, ' nnz=', nnz

  ! Run hungarian_match on the pre-built cost graph
  call hungarian_match(n, n, ptr, row_idx, cost, match_perm, matched, dualu, &
       dualv, stat)

  write(0, '(a,i8)') 'matched = ', matched
  write(0, '(a,i5)') 'stat = ', stat

  ! Check dual feasibility: u[i] + v[j] <= c[i,j] + eps
  max_infeasibility = 0.0_wp
  do j = 1, n
    do k = ptr(j), ptr(j+1)-1
      i = row_idx(k)
      v = dualu(i) + dualv(j) - cost(k)
      if (v > max_infeasibility) max_infeasibility = v
    end do
  end do
  write(0, '(a,es12.4)') 'max dual infeasibility (u+v-c): ', max_infeasibility

  ! Dump first few duals for comparison
  write(0, '(a)') 'First 10 dualu:'
  do i = 1, min(10, n)
    write(0, '(a,i5,a,es24.16)') '  u[', i, '] = ', dualu(i)
  end do
  write(0, '(a)') 'First 10 dualv:'
  do i = 1, min(10, n)
    write(0, '(a,i5,a,es24.16)') '  v[', i, '] = ', dualv(i)
  end do

end program spral_match_costgraph
