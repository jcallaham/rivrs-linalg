! Benchmark driver for SPRAL's SSIDS with internal ordering (match_order_metis).
!
! Unlike spral_full_solve.f90 which accepts a user-supplied ordering,
! this driver lets SPRAL compute its own MC64 matching + METIS ordering
! internally. This enables fair timing comparisons where each solver
! uses its own ordering pipeline.
!
! Input (stdin or file):
!   n nnz_lower
!   col_ptr[1] .. col_ptr[n+1]    (1-indexed, one per line)
!   row[1] val[1]                  (one per line, row 1-indexed)
!   ...
!   row[nnz_lower] val[nnz_lower]
!
! Output (stdout):
!   Structured key-value pairs between SPRAL_BENCHMARK_BEGIN / END sentinels.
!   All diagnostics go to stderr.
!
! OpenMP thread count is controlled externally via OMP_NUM_THREADS.
!
! Compile (after running tools/build_spral.sh):
!   gfortran -O2 -fopenmp -I /tmp/spral_ssids -o /tmp/spral_benchmark \
!     tools/spral_benchmark.f90 \
!     -Wl,--whole-archive /tmp/spral_ssids/libspral.a -Wl,--no-whole-archive \
!     $METIS_LIB -lopenblas -lstdc++ -lm -lgomp

program spral_benchmark_driver
  use spral_ssids
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, k, iunit, nargs, ios
  integer, allocatable :: ptr(:), row(:)
  real(wp), allocatable :: val(:), rhs(:), x(:)
  character(len=256) :: filename

  type(ssids_akeep) :: akeep
  type(ssids_fkeep) :: fkeep
  type(ssids_options) :: options
  type(ssids_inform) :: inform

  ! Timing
  integer(kind=8) :: t_start, t_end, count_rate
  real(wp) :: analyse_s, factor_s, solve_s

  ! Backward error computation
  real(wp), allocatable :: ax(:)
  real(wp) :: norm_a, norm_x, norm_b, norm_r, be
  integer :: j
  integer(kind=8) :: col_start, col_end

  ! Determine input source: file argument or stdin
  nargs = command_argument_count()
  if (nargs >= 1) then
    call get_command_argument(1, filename)
    iunit = 20
    open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(0, '(a,a)') 'ERROR: Cannot open file: ', trim(filename)
      stop 1
    end if
    write(0, '(a,a)') 'Reading matrix from file: ', trim(filename)
  else
    iunit = 5  ! stdin
    write(0, '(a)') 'Reading matrix from stdin'
  end if

  call system_clock(count_rate=count_rate)

  ! Read dimensions
  read(iunit, *) n, nnz

  allocate(ptr(n+1), row(nnz), val(nnz))
  allocate(rhs(n), x(n), ax(n))

  ! Read column pointers (1-indexed, lower triangle)
  do i = 1, n+1
    read(iunit, *) ptr(i)
  end do

  ! Read row indices and values
  do k = 1, nnz
    read(iunit, *) row(k), val(k)
  end do

  if (iunit /= 5) close(iunit)

  write(0, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Generate RHS: b = A * ones(n)
  ! A is stored as lower triangle; mirror for symmetric matvec
  rhs(1:n) = 0.0_wp
  do j = 1, n
    col_start = ptr(j)
    col_end = ptr(j+1) - 1
    do k = int(col_start), int(col_end)
      i = row(k)
      ! Lower triangle: a(i,j)
      rhs(i) = rhs(i) + val(k)
      if (i /= j) then
        ! Upper triangle (symmetric): a(j,i)
        rhs(j) = rhs(j) + val(k)
      end if
    end do
  end do

  ! Configure SPRAL SSIDS — let it do its own ordering
  options%ordering = 2        ! match_order_metis (MC64 + METIS)
  options%scaling = 3         ! use scaling from match_order_metis
  options%action = .true.     ! continue on singularity
  options%print_level = 0     ! quiet (diagnostics to unit 6 suppressed)
  options%ignore_numa = .false.
  options%u = 0.01_wp         ! default threshold parameter

  ! Analyse phase (SPRAL computes ordering internally)
  call system_clock(t_start)
  call ssids_analyse(.true., n, ptr, row, akeep, options, inform, val=val)
  call system_clock(t_end)
  analyse_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (inform%flag < 0) then
    write(0, '(a,i5)') 'ssids_analyse error flag=', inform%flag
    write(*, '(a)') 'SPRAL_BENCHMARK_BEGIN'
    write(*, '(a,i8)') 'n ', n
    write(*, '(a,i10)') 'nnz_lower ', nnz
    write(*, '(a,i5)') 'analyse_flag ', inform%flag
    write(*, '(a,i5)') 'factor_flag -999'
    write(*, '(a,i5)') 'solve_flag -999'
    write(*, '(a)') 'SPRAL_BENCHMARK_END'
    stop 1
  end if
  write(0, '(a,i5,a,f8.3,a)') 'analyse flag=', inform%flag, ' (', analyse_s, 's)'

  ! Factor phase (indefinite: posdef=.false.)
  call system_clock(t_start)
  call ssids_factor(.false., val, akeep, fkeep, options, inform)
  call system_clock(t_end)
  factor_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (inform%flag < 0) then
    write(0, '(a,i5)') 'ssids_factor error flag=', inform%flag
    write(*, '(a)') 'SPRAL_BENCHMARK_BEGIN'
    write(*, '(a,i8)') 'n ', n
    write(*, '(a,i10)') 'nnz_lower ', nnz
    write(*, '(a,i5)') 'analyse_flag ', 0
    write(*, '(a,i5)') 'factor_flag ', inform%flag
    write(*, '(a,i5)') 'solve_flag -999'
    write(*, '(a)') 'SPRAL_BENCHMARK_END'
    stop 1
  end if
  write(0, '(a,i5,a,f8.3,a)') 'factor flag=', inform%flag, ' (', factor_s, 's)'

  ! Solve phase
  x(1:n) = rhs(1:n)
  call system_clock(t_start)
  call ssids_solve(x, akeep, fkeep, options, inform)
  call system_clock(t_end)
  solve_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (inform%flag < 0) then
    write(0, '(a,i5)') 'ssids_solve error flag=', inform%flag
  end if
  write(0, '(a,i5,a,f8.3,a)') 'solve flag=', inform%flag, ' (', solve_s, 's)'

  ! Compute backward error: ||Ax - b||_2 / (||A||_F * ||x||_2 + ||b||_2)
  ! Uses L2/Frobenius norms to match rivrs sparse_backward_error() and
  ! the convention in Duff, Hogg & Lopez (2020).
  !
  ! A is stored as lower triangle CSC, but represents full symmetric matrix.
  ! For ||A||_F: lower-triangle storage means each off-diagonal entry appears
  ! once, but in the full matrix it appears twice. So:
  !   ||A||_F^2 = sum(diag^2) + 2*sum(off_diag^2)
  ! However, rivrs stores FULL symmetric CSC (both triangles), so its
  ! sparse_backward_error sums all stored entries: sum(v^2) which equals
  !   sum(diag^2) + 2*sum(off_diag^2) = ||A||_F^2
  ! We replicate that here from the lower triangle.
  ax(1:n) = 0.0_wp
  norm_a = 0.0_wp
  do j = 1, n
    col_start = ptr(j)
    col_end = ptr(j+1) - 1
    do k = int(col_start), int(col_end)
      i = row(k)
      ! Lower triangle: a(i,j)
      ax(i) = ax(i) + val(k) * x(j)
      ! Frobenius norm: diagonal counted once, off-diagonal counted twice
      if (i == j) then
        norm_a = norm_a + val(k)**2
      else
        norm_a = norm_a + 2.0_wp * val(k)**2
        ! Upper triangle (symmetric): a(j,i)
        ax(j) = ax(j) + val(k) * x(i)
      end if
    end do
  end do
  norm_a = sqrt(norm_a)

  ! Residual: ||Ax - b||_2, ||x||_2, ||b||_2
  norm_r = 0.0_wp
  norm_x = 0.0_wp
  norm_b = 0.0_wp
  do i = 1, n
    norm_r = norm_r + (ax(i) - rhs(i))**2
    norm_x = norm_x + x(i)**2
    norm_b = norm_b + rhs(i)**2
  end do
  norm_r = sqrt(norm_r)
  norm_x = sqrt(norm_x)
  norm_b = sqrt(norm_b)

  if (norm_a * norm_x + norm_b > 0.0_wp) then
    be = norm_r / (norm_a * norm_x + norm_b)
  else
    be = 0.0_wp
  end if

  write(0, '(a,es12.4)') 'Backward error: ', be
  write(0, '(a,i8,a,i8,a,i8,a,i8)') &
    'delays=', inform%num_delay, &
    ' 2x2=', inform%num_two, &
    ' not_1st=', inform%not_first_pass, &
    ' not_2nd=', inform%not_second_pass

  ! Output structured results between sentinels
  write(*, '(a)') 'SPRAL_BENCHMARK_BEGIN'
  write(*, '(a,i8)') 'n ', n
  write(*, '(a,i10)') 'nnz_lower ', nnz
  write(*, '(a,i5)') 'analyse_flag ', 0
  write(*, '(a,i5)') 'factor_flag ', inform%flag
  write(*, '(a,i5)') 'solve_flag ', 0
  write(*, '(a,es24.16)') 'analyse_s ', analyse_s
  write(*, '(a,es24.16)') 'factor_s ', factor_s
  write(*, '(a,es24.16)') 'solve_s ', solve_s
  write(*, '(a,es24.16)') 'backward_error ', be
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
  write(*, '(a,i8)') 'matrix_rank ', inform%matrix_rank
  write(*, '(a)') 'SPRAL_BENCHMARK_END'

  ! Cleanup
  call ssids_free(akeep, fkeep, i)

end program spral_benchmark_driver
