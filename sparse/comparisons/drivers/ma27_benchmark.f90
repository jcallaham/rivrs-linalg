! Benchmark driver for HSL MA27 on symmetric indefinite matrices.
!
! MA27 is the classic multifrontal solver by Duff & Reid (1982).
! It uses its own built-in minimum degree ordering (not METIS).
!
! Follows the same subprocess protocol as spral_benchmark.f90:
! reads a matrix, runs analyse/factor/solve, outputs sentinel-delimited
! key-value results for parsing by the Rust orchestration binary.
!
! Input (stdin or file):
!   n nnz_lower
!   row[1] col[1] val[1]     (COO, 1-indexed, lower triangle)
!   ...
!   row[nnz] col[nnz] val[nnz]
!
! Output (stdout):
!   Structured key-value pairs between MA27_BENCHMARK_BEGIN / END sentinels.
!   All diagnostics go to stderr.
!
! Compile (requires MA27 source from HSL):
!   comparisons/drivers/build_ma27.sh

program ma27_benchmark_driver
  implicit none

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, k, iunit, nargs, ios, j
  integer, allocatable :: irn(:), icn(:)
  real(wp), allocatable :: avals(:)
  character(len=256) :: filename

  ! MA27 workspace
  integer :: liw, la, iflag, nsteps, maxfrt
  integer, allocatable :: iw(:), ikeep(:), iw1(:)
  real(wp), allocatable :: a_factor(:), w(:), rhs(:), x(:)
  integer :: icntl(30), info(20)
  real(wp) :: cntl(5), ops

  ! Timing
  integer(kind=8) :: t_start, t_end, count_rate
  real(wp) :: analyse_s, factor_s, solve_s

  ! Backward error computation
  real(wp), allocatable :: ax(:)
  real(wp) :: norm_a, norm_x, norm_b, norm_r, be

  ! Determine input source: file argument or stdin
  nargs = command_argument_count()
  if (nargs >= 1) then
    call get_command_argument(1, filename)
    iunit = 20
    open(unit=iunit, file=trim(filename), status='old', action='read', &
         iostat=ios)
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

  allocate(irn(nnz), icn(nnz), avals(nnz))
  allocate(rhs(n), x(n), ax(n))

  ! Read COO entries (1-indexed, lower triangle)
  do k = 1, nnz
    read(iunit, *) irn(k), icn(k), avals(k)
  end do

  if (iunit /= 5) close(iunit)

  write(0, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Generate RHS: b = A * ones(n)
  rhs(1:n) = 0.0_wp
  do k = 1, nnz
    i = irn(k)
    j = icn(k)
    rhs(i) = rhs(i) + avals(k)
    if (i /= j) then
      rhs(j) = rhs(j) + avals(k)
    end if
  end do

  ! Initialize MA27 control parameters
  call ma27id(icntl, cntl)
  ! Suppress output
  icntl(1) = 0  ! error message stream (0 = suppress)
  icntl(2) = 0  ! diagnostic message stream

  ! Allocate workspace for analysis
  ! MA27 documentation recommends liw >= 2*nnz + 3*n + 1
  liw = 3 * nnz + 3 * n + 1
  allocate(iw(liw))
  allocate(ikeep(3*n))
  allocate(iw1(2*n))

  ! Analysis phase (MA27AD)
  call system_clock(t_start)
  call ma27ad(n, nnz, irn, icn, iw, liw, ikeep, iw1, nsteps, iflag, &
              icntl, cntl, info, ops)
  call system_clock(t_end)
  analyse_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (iflag /= 0) then
    write(0, '(a,i5)') 'MA27AD warning/error iflag=', iflag
  end if
  if (info(1) < 0) then
    write(0, '(a,i5)') 'MA27AD error info(1)=', info(1)
    write(*, '(a)') 'MA27_BENCHMARK_BEGIN'
    write(*, '(a,i5)') 'analyse_flag ', info(1)
    write(*, '(a)') 'factor_flag -999'
    write(*, '(a)') 'solve_flag -999'
    write(*, '(a)') 'MA27_BENCHMARK_END'
    stop 1
  end if
  write(0, '(a,i5,a,f8.3,a)') 'analyse iflag=', iflag, ' (', analyse_s, 's)'

  ! Factorization phase (MA27BD) with retry loop on workspace-too-small.
  ! MA27BD returns info(1) = -3 (la too small) or -4 (liw too small).
  ! On failure, info(5)/info(6) report the required sizes. We reallocate
  ! with a safety margin and retry, up to 3 attempts.

  ! Initial workspace allocation from analysis estimates
  la = max(int(info(5) * 2.0), nnz)
  deallocate(iw)
  liw = max(int(info(6) * 2.0), nnz)
  allocate(iw(liw))
  allocate(a_factor(la))

  do i = 1, 3
    ! Copy matrix values into factor workspace (MA27BD overwrites a_factor)
    a_factor(1:nnz) = avals(1:nnz)

    ! Re-run analysis if iw was reallocated (MA27 requires this)
    if (i > 1) then
      call ma27ad(n, nnz, irn, icn, iw, liw, ikeep, iw1, nsteps, iflag, &
                  icntl, cntl, info, ops)
      ! Grow la if analysis now reports larger estimate
      if (la < info(5)) then
        la = int(info(5) * 3.0)
        deallocate(a_factor)
        allocate(a_factor(la))
        a_factor(1:nnz) = avals(1:nnz)
      end if
    end if

    call system_clock(t_start)
    call ma27bd(n, nnz, irn, icn, a_factor, la, iw, liw, ikeep, nsteps, &
                maxfrt, iw1, icntl, cntl, info)
    call system_clock(t_end)

    if (info(1) /= -3 .and. info(1) /= -4) exit  ! success or fatal error

    write(0, '(a,i5,a,i1)') 'MA27BD workspace too small (info(1)=', &
      info(1), '), retry ', i

    ! Grow the workspace that was too small.
    ! MA27BD: -3 = LIW too small (INFO(2) = required LIW)
    !         -4 = LA too small  (INFO(2) = required LA)
    if (info(1) == -3) then
      liw = max(int(info(2) * 2.0), liw * 3)
      deallocate(iw)
      allocate(iw(liw))
    end if
    if (info(1) == -4) then
      la = max(int(info(2) * 2.0), la * 3)
      deallocate(a_factor)
      allocate(a_factor(la))
    end if
  end do
  factor_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (info(1) < 0) then
    write(0, '(a,i5)') 'MA27BD error info(1)=', info(1)
    write(*, '(a)') 'MA27_BENCHMARK_BEGIN'
    write(*, '(a,es24.16)') 'analyse_s ', analyse_s
    write(*, '(a,i5)') 'analyse_flag ', 0
    write(*, '(a,i5)') 'factor_flag ', info(1)
    write(*, '(a)') 'solve_flag -999'
    write(*, '(a)') 'MA27_BENCHMARK_END'
    stop 1
  end if
  write(0, '(a,i5,a,f8.3,a)') 'factor info(1)=', info(1), &
    ' (', factor_s, 's)'

  ! Solve phase (MA27CD)
  allocate(w(maxfrt))
  x(1:n) = rhs(1:n)

  call system_clock(t_start)
  call ma27cd(n, a_factor, la, iw, liw, w, maxfrt, x, iw1, nsteps, &
              icntl, info)
  call system_clock(t_end)
  solve_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (info(1) < 0) then
    write(0, '(a,i5)') 'MA27CD error info(1)=', info(1)
  end if
  write(0, '(a,i5,a,f8.3,a)') 'solve info(1)=', info(1), &
    ' (', solve_s, 's)'

  ! Compute backward error: ||Ax - b||_2 / (||A||_F * ||x||_2 + ||b||_2)
  ax(1:n) = 0.0_wp
  norm_a = 0.0_wp
  do k = 1, nnz
    i = irn(k)
    j = icn(k)
    ax(i) = ax(i) + avals(k) * x(j)
    if (i == j) then
      norm_a = norm_a + avals(k)**2
    else
      norm_a = norm_a + 2.0_wp * avals(k)**2
      ax(j) = ax(j) + avals(k) * x(i)
    end if
  end do
  norm_a = sqrt(norm_a)

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

  ! Output structured results between sentinels
  write(*, '(a)') 'MA27_BENCHMARK_BEGIN'
  write(*, '(a,i8)') 'n ', n
  write(*, '(a,i10)') 'nnz_lower ', nnz
  write(*, '(a,i5)') 'analyse_flag ', 0
  write(*, '(a,i5)') 'factor_flag ', 0
  write(*, '(a,i5)') 'solve_flag ', info(1)
  write(*, '(a,es24.16)') 'analyse_s ', analyse_s
  write(*, '(a,es24.16)') 'factor_s ', factor_s
  write(*, '(a,es24.16)') 'solve_s ', solve_s
  write(*, '(a,es24.16)') 'backward_error ', be
  write(*, '(a,i8)') 'maxfrt ', maxfrt
  write(*, '(a,i8)') 'nsteps ', nsteps
  write(*, '(a,es24.16)') 'est_ops ', ops
  write(*, '(a,i10)') 'factor_size ', info(9)
  write(*, '(a,i8)') 'num_2x2 ', info(14)
  write(*, '(a,i8)') 'num_neg ', info(15)
  write(*, '(a)') 'MA27_BENCHMARK_END'

end program ma27_benchmark_driver
