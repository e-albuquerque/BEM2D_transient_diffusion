!-----------------------------------------------------------------------
! MODULE: solver_mod
!> @brief Linear system solver module using LAPACK (DGESV, DGETRF, DGETRS).
!-----------------------------------------------------------------------
module solver_mod
    implicit none

    ! Explicit interfaces for LAPACK double-precision routines
    interface
        subroutine DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            implicit none
            integer, intent(in)    :: N, NRHS, LDA, LDB
            integer, intent(out)   :: IPIV(N), INFO
            real,    intent(inout) :: A(LDA, *), B(LDB, *)
        end subroutine DGESV

        subroutine DGETRF(M, N, A, LDA, IPIV, INFO)
            implicit none
            integer, intent(in)    :: M, N, LDA
            integer, intent(out)   :: IPIV(*), INFO
            real,    intent(inout) :: A(LDA, *)
        end subroutine DGETRF

        subroutine DGETRS(TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO)
            implicit none
            character(len=1), intent(in) :: TRANS
            integer, intent(in)    :: N, NRHS, LDA, LDB
            integer, intent(in)    :: IPIV(*)
            real,    intent(in)    :: A(LDA, *)
            real,    intent(inout) :: B(LDB, *)
            integer, intent(out)   :: INFO
        end subroutine DGETRS
    end interface

contains

    !-------------------------------------------------------------------
    !> @brief Computes LU factorization of system matrix A using LAPACK
    !>        DGETRF so it can be reused across all time steps.
    !-------------------------------------------------------------------
    subroutine FACTOR_SYSTEM(N, A, LDA, IPIV, IPR)
        implicit none
        integer, intent(in)    :: N, LDA, IPR
        real,    intent(inout) :: A(LDA, LDA)
        integer, intent(out)   :: IPIV(N)
        integer :: INFO

        call DGETRF(N, N, A, LDA, IPIV, INFO)
        if (INFO /= 0) then
            write(IPR, 10) INFO
            write(*,   10) INFO
            stop 1
        end if
10      format(' **** SINGULARITY OR LAPACK DGETRF ERROR (INFO =', I5, ')')
    end subroutine FACTOR_SYSTEM

    !-------------------------------------------------------------------
    !> @brief Solves A * X = B using precomputed LU factors from DGETRF
    !>        via LAPACK DGETRS.
    !-------------------------------------------------------------------
    subroutine SOLVE_FACTORED(N, A, LDA, IPIV, B, IPR)
        implicit none
        integer, intent(in)    :: N, LDA, IPR
        real,    intent(in)    :: A(LDA, LDA)
        integer, intent(in)    :: IPIV(N)
        real,    intent(inout) :: B(LDA)
        integer :: INFO

        call DGETRS('N', N, 1, A, LDA, IPIV, B, LDA, INFO)
        if (INFO /= 0) then
            write(IPR, 20) INFO
            write(*,   20) INFO
            stop 1
        end if
20      format(' **** LAPACK DGETRS ERROR (INFO =', I5, ')')
    end subroutine SOLVE_FACTORED

    !-------------------------------------------------------------------
    !> @brief Direct one-shot solver A * X = B using LAPACK DGESV
    !>        (drop-in replacement for legacy SLNPD).
    !-------------------------------------------------------------------
    subroutine SLNPD(MX, NN, A, B, IPR)
        implicit none
        integer, intent(in)    :: MX, NN, IPR
        real,    intent(inout) :: A(MX, MX), B(MX)
        integer, allocatable   :: IPIV_LOC(:)
        integer :: INFO

        allocate(IPIV_LOC(NN))
        call DGESV(NN, 1, A, MX, IPIV_LOC, B, MX, INFO)
        if (INFO /= 0) then
            write(IPR, 30) INFO
            write(*,   30) INFO
        end if
        deallocate(IPIV_LOC)
30      format(' **** SINGULARITY OR LAPACK DGESV ERROR (INFO =', I5, ')')
    end subroutine SLNPD

end module solver_mod
