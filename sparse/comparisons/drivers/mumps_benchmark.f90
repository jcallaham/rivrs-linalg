! Benchmark driver for MUMPS on symmetric indefinite matrices.
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
!   Structured key-value pairs between MUMPS_BENCHMARK_BEGIN / END sentinels.
!   All diagnostics go to stderr.
!
! Environment variables:
!   MUMPS_ORDERING: auto (default), metis, amd, scotch, pord
!     Maps to ICNTL(7) values: 7, 5, 0, 3, 4
!   OPENBLAS_NUM_THREADS / OMP_NUM_THREADS: control BLAS threading
!     (set to 1 for single-threaded benchmarking)
!
! Compile (after ensuring libmumps-seq-dev is installed):
!   comparisons/drivers/build_mumps.sh

program mumps_benchmark_driver
  implicit none
  include 'dmumps_struc.h'

  integer, parameter :: wp = kind(0d0)
  integer :: n, nnz, i, k, iunit, nargs, ios
  integer, allocatable, target :: irn(:), jcn(:)
  real(wp), allocatable, target :: a_vals(:), rhs(:)
  real(wp), allocatable :: x(:)
  character(len=256) :: filename, ordering_str

  type(dmumps_struc) :: mumps

  ! Timing
  integer(kind=8) :: t_start, t_end, count_rate
  real(wp) :: analyse_s, factor_s, solve_s

  ! Backward error computation
  real(wp), allocatable :: ax(:)
  real(wp) :: norm_a, norm_x, norm_b, norm_r, be
  integer :: j

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

  allocate(irn(nnz), jcn(nnz), a_vals(nnz))
  allocate(rhs(n), x(n), ax(n))

  ! Read COO entries (1-indexed, lower triangle)
  do k = 1, nnz
    read(iunit, *) irn(k), jcn(k), a_vals(k)
  end do

  if (iunit /= 5) close(iunit)

  write(0, '(a,i8,a,i10)') 'Matrix: n=', n, ' nnz_lower=', nnz

  ! Generate RHS: b = A * ones(n)
  ! A is stored as lower triangle COO; mirror for symmetric matvec
  rhs(1:n) = 0.0_wp
  do k = 1, nnz
    i = irn(k)
    j = jcn(k)
    ! Lower triangle: a(i,j)
    rhs(i) = rhs(i) + a_vals(k)
    if (i /= j) then
      ! Upper triangle (symmetric): a(j,i)
      rhs(j) = rhs(j) + a_vals(k)
    end if
  end do

  ! Initialize MUMPS
  mumps%COMM = 0  ! not used for sequential (dummy MPI)
  mumps%SYM = 2   ! symmetric indefinite (LDL^T with 2x2 pivots)
  mumps%PAR = 1   ! host participates in computation
  mumps%JOB = -1  ! initialize

  call dmumps(mumps)

  if (mumps%INFOG(1) < 0) then
    write(0, '(a,i5)') 'MUMPS init error INFOG(1)=', mumps%INFOG(1)
    write(*, '(a)') 'MUMPS_BENCHMARK_BEGIN'
    write(*, '(a,i5)') 'error ', mumps%INFOG(1)
    write(*, '(a)') 'MUMPS_BENCHMARK_END'
    stop 1
  end if

  ! Configure MUMPS
  ! Suppress stdout output
  mumps%ICNTL(1) = 0  ! error messages stream (0 = suppress)
  mumps%ICNTL(2) = 0  ! diagnostic messages stream
  mumps%ICNTL(3) = 0  ! global info stream
  mumps%ICNTL(4) = 0  ! print level (0 = no output)

  ! Memory relaxation for hard matrices
  mumps%ICNTL(14) = 50  ! 50% extra workspace

  ! Ordering selection
  mumps%ICNTL(7) = 7  ! default: auto (MUMPS picks best)
  call get_environment_variable('MUMPS_ORDERING', ordering_str, status=ios)
  if (ios == 0) then
    select case (trim(ordering_str))
      case ('auto')
        mumps%ICNTL(7) = 7
      case ('amd')
        mumps%ICNTL(7) = 0
      case ('scotch')
        mumps%ICNTL(7) = 3
      case ('pord')
        mumps%ICNTL(7) = 4
      case ('metis')
        mumps%ICNTL(7) = 5
      case default
        write(0, '(a,a)') 'Warning: unknown MUMPS_ORDERING=', &
          trim(ordering_str)
    end select
  end if

  ! Set matrix data (pointer assignment — arrays have TARGET attribute)
  mumps%N = n
  mumps%NZ = nnz
  mumps%IRN => irn
  mumps%JCN => jcn
  mumps%A => a_vals
  mumps%RHS => rhs

  ! Analysis phase
  mumps%JOB = 1
  call system_clock(t_start)
  call dmumps(mumps)
  call system_clock(t_end)
  analyse_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (mumps%INFOG(1) < 0) then
    write(0, '(a,i5)') 'MUMPS analyse error INFOG(1)=', mumps%INFOG(1)
    write(*, '(a)') 'MUMPS_BENCHMARK_BEGIN'
    write(*, '(a,i5)') 'analyse_flag ', mumps%INFOG(1)
    write(*, '(a)') 'factor_flag -999'
    write(*, '(a)') 'solve_flag -999'
    write(*, '(a)') 'MUMPS_BENCHMARK_END'
    ! Cleanup
    mumps%JOB = -2
    call dmumps(mumps)
    stop 1
  end if
  write(0, '(a,i5,a,f8.3,a)') 'analyse INFOG(1)=', mumps%INFOG(1), &
    ' (', analyse_s, 's)'

  ! Factorization phase
  mumps%JOB = 2
  call system_clock(t_start)
  call dmumps(mumps)
  call system_clock(t_end)
  factor_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (mumps%INFOG(1) < 0) then
    write(0, '(a,i5)') 'MUMPS factor error INFOG(1)=', mumps%INFOG(1)
    write(*, '(a)') 'MUMPS_BENCHMARK_BEGIN'
    write(*, '(a,es24.16)') 'analyse_s ', analyse_s
    write(*, '(a,i5)') 'analyse_flag ', 0
    write(*, '(a,i5)') 'factor_flag ', mumps%INFOG(1)
    write(*, '(a)') 'solve_flag -999'
    write(*, '(a)') 'MUMPS_BENCHMARK_END'
    mumps%JOB = -2
    call dmumps(mumps)
    stop 1
  end if
  write(0, '(a,i5,a,f8.3,a)') 'factor INFOG(1)=', mumps%INFOG(1), &
    ' (', factor_s, 's)'

  ! Solve phase (RHS is overwritten with solution in-place)
  mumps%JOB = 3
  call system_clock(t_start)
  call dmumps(mumps)
  call system_clock(t_end)
  solve_s = real(t_end - t_start, wp) / real(count_rate, wp)

  if (mumps%INFOG(1) < 0) then
    write(0, '(a,i5)') 'MUMPS solve error INFOG(1)=', mumps%INFOG(1)
  end if
  write(0, '(a,i5,a,f8.3,a)') 'solve INFOG(1)=', mumps%INFOG(1), &
    ' (', solve_s, 's)'

  ! Solution is in mumps%RHS (which points to rhs array)
  x(1:n) = rhs(1:n)

  ! Regenerate RHS for backward error (MUMPS overwrote it)
  rhs(1:n) = 0.0_wp
  do k = 1, nnz
    i = irn(k)
    j = jcn(k)
    rhs(i) = rhs(i) + a_vals(k)
    if (i /= j) then
      rhs(j) = rhs(j) + a_vals(k)
    end if
  end do

  ! Compute backward error: ||Ax - b||_2 / (||A||_F * ||x||_2 + ||b||_2)
  ax(1:n) = 0.0_wp
  norm_a = 0.0_wp
  do k = 1, nnz
    i = irn(k)
    j = jcn(k)
    ! Lower triangle: a(i,j)
    ax(i) = ax(i) + a_vals(k) * x(j)
    if (i == j) then
      norm_a = norm_a + a_vals(k)**2
    else
      norm_a = norm_a + 2.0_wp * a_vals(k)**2
      ! Upper triangle (symmetric): a(j,i)
      ax(j) = ax(j) + a_vals(k) * x(i)
    end if
  end do
  norm_a = sqrt(norm_a)

  ! Residual norms
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
  write(*, '(a)') 'MUMPS_BENCHMARK_BEGIN'
  write(*, '(a,i8)') 'n ', n
  write(*, '(a,i10)') 'nnz_lower ', nnz
  write(*, '(a,i5)') 'analyse_flag ', 0
  write(*, '(a,i5)') 'factor_flag ', 0
  write(*, '(a,i5)') 'solve_flag ', mumps%INFOG(1)
  write(*, '(a,es24.16)') 'analyse_s ', analyse_s
  write(*, '(a,es24.16)') 'factor_s ', factor_s
  write(*, '(a,es24.16)') 'solve_s ', solve_s
  write(*, '(a,es24.16)') 'backward_error ', be
  write(*, '(a,i8)') 'num_neg ', mumps%INFOG(12)
  write(*, '(a,i8)') 'ordering_used ', mumps%INFOG(7)
  write(*, '(a,i20)') 'factor_entries ', mumps%INFOG(29)
  write(*, '(a,es24.16)') 'est_flops ', mumps%RINFOG(3)
  write(*, '(a)') 'MUMPS_BENCHMARK_END'

  ! Cleanup
  mumps%JOB = -2
  call dmumps(mumps)

end program mumps_benchmark_driver
