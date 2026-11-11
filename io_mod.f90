!-----------------------------------------------------------------------
! MODULE: io_mod
!> @brief Input and output routines for modernized 2D Transient Diffusion
!>        BEM (DCC).
!-----------------------------------------------------------------------
module io_mod
    use data_mod
    implicit none

    character(len=80) :: JOB_TITLE

contains

    !-------------------------------------------------------------------
    !> @brief Reads problem parameters, mesh nodes, element connectivity,
    !>        time-varying BCs, and internal points from dcc_modern.dat,
    !>        and dynamically allocates arrays.
    !-------------------------------------------------------------------
    subroutine INPUT_DIFF()
        implicit none
        integer :: I, KN, KE, IT
        real, allocatable :: VAL_TMP(:)

        read(INP, '(A)') JOB_TITLE
        read(INP, *) NE, NN, NINTER, NPI, K_DIFF, DELTA_T, NG

        if (NG < 2) NG = 8

        call allocate_arrays()

        ! Read geometric boundary node coordinates (1..NN)
        do I = 1, NN
            read(INP, *) KN, XG(KN), YG(KN)
        end do

        ! Read element connectivity (No1, No2) and BC code (0=U, 1=P)
        do I = 1, NE
            read(INP, *) KE, CONEC(KE, 1), CONEC(KE, 2), KODE(KE)
        end do

        ! Read prescribed boundary values for each element across NINTER time steps
        allocate(VAL_TMP(NINTER))
        do I = 1, NE
            read(INP, *) (VAL_TMP(IT), IT = 1, NINTER)
            if (KODE(I) == 0) then
                FI(I, 1:NINTER) = VAL_TMP(1:NINTER)
            else
                DFI(I, 1:NINTER) = VAL_TMP(1:NINTER)
            end if
        end do
        deallocate(VAL_TMP)

        ! Read internal point coordinates
        if (NPI > 0) then
            do I = 1, NPI
                read(INP, *) KN, XI(KN), YI(KN)
            end do
        end if
    end subroutine INPUT_DIFF

    !-------------------------------------------------------------------
    !> @brief Exports formatted report (dcc_modern.out) and tabular data
    !>        files (ubou.dat, pbou.dat, uint.dat) for Python postprocessing.
    !-------------------------------------------------------------------
    subroutine OUTPT_DIFF(ELAPSED_SEC, NTHREADS)
        implicit none
        real,    intent(in) :: ELAPSED_SEC
        integer, intent(in) :: NTHREADS
        integer :: I, IT, KN

        write(IPR, '(A)') repeat('=', 78)
        write(IPR, '(A)') trim(JOB_TITLE)
        write(IPR, '(A)') repeat('=', 78)
        write(IPR, '(A, I6)')     ' Boundary Elements (NE)    : ', NE
        write(IPR, '(A, I6)')     ' Geometric Nodes (NN)      : ', NN
        write(IPR, '(A, I6)')     ' Time Steps (NINTER)       : ', NINTER
        write(IPR, '(A, I6)')     ' Internal Points (NPI)     : ', NPI
        write(IPR, '(A, I6)')     ' Gauss Points (NG)         : ', NG
        write(IPR, '(A, E14.6)')  ' Diffusivity (K)           : ', K_DIFF
        write(IPR, '(A, E14.6)')  ' Time Step (DeltaT)        : ', DELTA_T
        write(IPR, '(A, I6)')     ' OpenMP Threads            : ', NTHREADS
        write(IPR, '(A, F10.4, A)') ' Wall Elapsed Time         : ', ELAPSED_SEC, ' s'
        write(IPR, '(A)') repeat('=', 78)

        do IT = 1, NINTER
            write(IPR, '(/, A, I5, A, F10.4)') ' Tempo: ', IT, '   t = ', TF(IT)
            write(IPR, '(A)') '  Elem.      Xm          Ym            U               P'
            do I = 1, NE
                write(IPR, '(1X, I5, 2(2X, F10.4), 2(2X, E14.6))') &
                    I, XM(I), YM(I), FI(I, IT), DFI(I, IT)
            end do

            if (NPI > 0) then
                write(IPR, '(A)') ' PInter.     Xi          Yi            U'
                do KN = 1, NPI
                    write(IPR, '(1X, I5, 2(2X, F10.4), 2X, E14.6)') &
                        KN, XI(KN), YI(KN), UPI(KN, IT)
                end do
            end if
        end do

        ! Export compact tables for Python postprocessing:
        ! ubou.dat: Xm Ym U(t1) ... U(t_NINTER)
        ! pbou.dat: Xm Ym P(t1) ... P(t_NINTER)
        ! uint.dat: Xi Yi U(t1) ... U(t_NINTER)
        do I = 1, NE
            write(IP1, '(*(ES22.14, 1X))') XM(I), YM(I), (FI(I, IT), IT = 1, NINTER)
            write(IP2, '(*(ES22.14, 1X))') XM(I), YM(I), (DFI(I, IT), IT = 1, NINTER)
        end do

        if (NPI > 0) then
            do KN = 1, NPI
                write(IP3, '(*(ES22.14, 1X))') XI(KN), YI(KN), (UPI(KN, IT), IT = 1, NINTER)
            end do
        end if
    end subroutine OUTPT_DIFF

end module io_mod
