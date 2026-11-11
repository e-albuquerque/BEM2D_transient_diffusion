!-----------------------------------------------------------------------
! PROGRAM: DCC_MODERN (Modernized Fortran 90 Version)
!> @brief Solves 2D transient heat diffusion equations using the
!>        Boundary Element Method with constant elements in space and
!>        time, OpenMP parallelization, dynamic allocation, and LAPACK.
!-----------------------------------------------------------------------
program DCC_MODERN
    use omp_lib
    use data_mod
    use io_mod
    use assembly_mod
    use solver_mod
    implicit none

    character(len=64) :: FILEIN, FILEOUT
    real :: T_START, T_END, ELAPSED
    integer :: NTHREADS

    ! Unit numbers
    INP = 10
    IPR = 11
    IP1 = 12
    IP2 = 13
    IP3 = 14

    FILEIN  = 'dcc_modern.dat'
    FILEOUT = 'dcc_modern.out'

    open(INP, file=trim(FILEIN),  status='OLD',     action='READ')
    open(IPR, file=trim(FILEOUT), status='REPLACE', action='WRITE')
    open(IP1, file='ubou.dat',    status='REPLACE', action='WRITE')
    open(IP2, file='pbou.dat',    status='REPLACE', action='WRITE')
    open(IP3, file='uint.dat',    status='REPLACE', action='WRITE')

    T_START  = omp_get_wtime()
    NTHREADS = omp_get_max_threads()

    ! 1. Read input data and allocate dynamic arrays
    call INPUT_DIFF()

    ! 2. Precompute element midpoints, lengths, outward normals, and times
    call PRECOMPUTE_GEOMETRY()

    ! 3. Compute G and H influence matrices in parallel via OpenMP
    call MATGH()

    ! 4. Factor system matrix once (DGETRF) and march in time (DGETRS)
    call SOLVE_TRANSIENT()

    ! 5. Compute potential at internal points in parallel via OpenMP
    call PTINT()

    T_END   = omp_get_wtime()
    ELAPSED = T_END - T_START

    ! 6. Export formatted report and data tables for Python postprocessing
    call OUTPT_DIFF(ELAPSED, NTHREADS)

    write(*, '(A)') trim(JOB_TITLE)
    write(*, '(A, I0, A, F8.4, A)') 'Completed with ', NTHREADS, &
        ' OpenMP threads in ', ELAPSED, ' s.'

    ! 7. Release dynamic memory and close files
    call deallocate_arrays()

    close(INP)
    close(IPR)
    close(IP1)
    close(IP2)
    close(IP3)

end program DCC_MODERN
