! Dump SPRAL's greedy matching state (after hungarian_init_heurisitic, before Dijkstra)
!
! Usage: ./spral_greedy_compare <cost_graph_file>
! Cost graph file format (1-indexed col_ptr, 1-indexed row):
!   Line 1: n nnz
!   Next n+1 lines: col_ptr values (1-indexed)
!   Next nnz lines: row_idx(1-indexed) cost

program spral_greedy_compare
   use spral_scaling, only : hungarian_init_heurisitic
   implicit none
   integer, parameter :: wp = kind(1.0d0)
   integer, parameter :: long = selected_int_kind(18)
   real(wp), parameter :: RINF = huge(0.0_wp)

   integer :: n, nnz, i, j, num
   integer(long) :: k
   integer(long), allocatable :: ptr(:), jperm(:), search_from(:), longwork(:)
   integer, allocatable :: row(:), iperm(:)
   real(wp), allocatable :: val(:), dualu(:), d(:)
   character(len=512) :: filename
   integer :: funit, unmatched_count

   ! Read arguments
   call get_command_argument(1, filename)
   if (len_trim(filename) == 0) then
      write(*,*) 'Usage: spral_greedy_compare <cost_graph_file>'
      stop 1
   end if

   ! Read cost graph (already 1-indexed from export)
   funit = 20
   open(unit=funit, file=trim(filename), status='old', action='read')
   read(funit, *) n, nnz
   allocate(ptr(n+1), row(nnz), val(nnz))
   do i = 1, n+1
      read(funit, *) ptr(i)
   end do
   do k = 1, nnz
      read(funit, *) row(k), val(k)
   end do
   close(funit)

   write(*,'(A,I0,A,I0)') 'n=', n, ' nnz=', nnz

   ! Run SPRAL's greedy initialization
   allocate(iperm(n), jperm(n), dualu(n), d(n), longwork(n), search_from(n))
   iperm = 0
   jperm = 0
   num = 0
   search_from = 0

   call hungarian_init_heurisitic(n, n, ptr, row, val, num, iperm, jperm, &
        dualu, d, longwork, search_from)

   write(*,'(A,I0,A,I0)') 'SPRAL greedy matched: ', num, ' / ', n

   ! Dump first 30 unmatched columns (0-indexed for comparison with our code)
   unmatched_count = 0
   write(*,'(A)') 'First 30 unmatched columns (0-indexed):'
   do j = 1, n
      if (jperm(j) == 0) then
         unmatched_count = unmatched_count + 1
         if (unmatched_count <= 30) write(*,'(A,I0)') '  col ', j-1
      end if
   end do
   write(*,'(A,I0)') 'Total unmatched columns: ', n - num

   ! Compare dual properties
   write(*,'(A,I0)') 'Rows with dualu > 0: ', count(dualu > 1e-15)
   write(*,'(A,I0)') 'Rows with dualu = 0: ', count(abs(dualu) <= 1e-15)
   write(*,'(A,I0)') 'Rows with dualu = RINF: ', count(dualu >= RINF)
   write(*,'(A,ES15.6)') 'min(dualu): ', minval(dualu)
   write(*,'(A,ES15.6)') 'max(dualu): ', maxval(dualu)

end program spral_greedy_compare
