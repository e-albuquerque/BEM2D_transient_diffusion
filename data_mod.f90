!-----------------------------------------------------------------------
! MODULE: data_mod
!> @brief Global variables and dynamic memory allocation for modernized
!>        2D Transient Diffusion BEM (DCC).
!-----------------------------------------------------------------------
module data_mod
    implicit none

    ! Mathematical constant
    real, parameter :: PI = 3.14159265358979323846

    ! Problem dimensions and file unit numbers
    integer :: NE, NN, NINTER, NPI, NG
    integer :: INP, IPR, IP1, IP2, IP3

    ! Physical and temporal parameters
    real :: K_DIFF, DELTA_T

    ! Dynamically allocated geometry, connectivity, and BCs
    real,    allocatable :: XG(:), YG(:)
    integer, allocatable :: CONEC(:,:)
    real,    allocatable :: XM(:), YM(:), EL_LEN(:), NX_VEC(:), NY_VEC(:)
    integer, allocatable :: KODE(:)
    real,    allocatable :: FI(:,:), DFI(:,:)
    real,    allocatable :: TF(:)

    ! Dynamically allocated BEM influence and system matrices
    real,    allocatable :: G(:,:,:), H(:,:,:)
    real,    allocatable :: A_MAT(:,:), B_MAT(:,:), C_MAT(:,:)
    real,    allocatable :: RHS(:), VPR(:)
    integer, allocatable :: IPIV(:)

    ! Dynamically allocated internal point coordinates and results
    real,    allocatable :: XI(:), YI(:)
    real,    allocatable :: UPI(:,:)

contains

    !-------------------------------------------------------------------
    !> @brief Allocates all dynamic arrays based on NE, NN, NINTER, NPI.
    !-------------------------------------------------------------------
    subroutine allocate_arrays()
        implicit none

        allocate(XG(NN), YG(NN))
        allocate(CONEC(NE, 2))
        allocate(XM(NE), YM(NE), EL_LEN(NE), NX_VEC(NE), NY_VEC(NE))
        allocate(KODE(NE))
        allocate(FI(NE, NINTER), DFI(NE, NINTER))
        allocate(TF(NINTER))

        allocate(G(NE, NE, NINTER), H(NE, NE, NINTER))
        allocate(A_MAT(NE, NE), B_MAT(NE, NE), C_MAT(NE, NE))
        allocate(RHS(NE), VPR(NE))
        allocate(IPIV(NE))

        if (NPI > 0) then
            allocate(XI(NPI), YI(NPI))
            allocate(UPI(NPI, NINTER))
        else
            allocate(XI(1), YI(1))
            allocate(UPI(1, 1))
        end if

        XG     = 0.0
        YG     = 0.0
        CONEC  = 0
        XM     = 0.0
        YM     = 0.0
        EL_LEN = 0.0
        NX_VEC = 0.0
        NY_VEC = 0.0
        KODE   = 0
        FI     = 0.0
        DFI    = 0.0
        TF     = 0.0
        G      = 0.0
        H      = 0.0
        A_MAT  = 0.0
        B_MAT  = 0.0
        C_MAT  = 0.0
        RHS    = 0.0
        VPR    = 0.0
        IPIV   = 0
        XI     = 0.0
        YI     = 0.0
        UPI    = 0.0
    end subroutine allocate_arrays

    !-------------------------------------------------------------------
    !> @brief Deallocates all dynamic arrays.
    !-------------------------------------------------------------------
    subroutine deallocate_arrays()
        implicit none

        if (allocated(XG))     deallocate(XG)
        if (allocated(YG))     deallocate(YG)
        if (allocated(CONEC))  deallocate(CONEC)
        if (allocated(XM))     deallocate(XM)
        if (allocated(YM))     deallocate(YM)
        if (allocated(EL_LEN)) deallocate(EL_LEN)
        if (allocated(NX_VEC)) deallocate(NX_VEC)
        if (allocated(NY_VEC)) deallocate(NY_VEC)
        if (allocated(KODE))   deallocate(KODE)
        if (allocated(FI))     deallocate(FI)
        if (allocated(DFI))    deallocate(DFI)
        if (allocated(TF))     deallocate(TF)
        if (allocated(G))      deallocate(G)
        if (allocated(H))      deallocate(H)
        if (allocated(A_MAT))  deallocate(A_MAT)
        if (allocated(B_MAT))  deallocate(B_MAT)
        if (allocated(C_MAT))  deallocate(C_MAT)
        if (allocated(RHS))    deallocate(RHS)
        if (allocated(VPR))    deallocate(VPR)
        if (allocated(IPIV))   deallocate(IPIV)
        if (allocated(XI))     deallocate(XI)
        if (allocated(YI))     deallocate(YI)
        if (allocated(UPI))    deallocate(UPI)
    end subroutine deallocate_arrays

end module data_mod
