!-----------------------------------------------------------------------
! MODULE: assembly_mod
!> @brief Parallel OpenMP assembly of BEM transient influence matrices,
!>        time-marching solver driver, and internal point evaluation.
!-----------------------------------------------------------------------
module assembly_mod
    use data_mod
    use solver_mod
    implicit none

contains

    !-------------------------------------------------------------------
    !> @brief Computes NPG-point Gauss-Legendre nodes QSI and weights W
    !>        on [-1, 1] in double precision.
    !-------------------------------------------------------------------
    subroutine GAUSS_LEGENDRE(NPG, QSI, W)
        implicit none
        integer, intent(in)  :: NPG
        real,    intent(out) :: QSI(NPG), W(NPG)

        real, parameter :: TOL = 1.0e-14
        integer :: I, J, M
        real :: Z, Z1, P1, P2, P3, PP

        M = (NPG + 1) / 2
        do I = 1, M
            Z = cos(PI * (real(I) - 0.25) / (real(NPG) + 0.5))
            do
                P1 = 1.0
                P2 = 0.0
                do J = 1, NPG
                    P3 = P2
                    P2 = P1
                    P1 = ((2.0 * real(J) - 1.0) * Z * P2 - (real(J) - 1.0) * P3) / real(J)
                end do
                PP = real(NPG) * (Z * P1 - P2) / (Z * Z - 1.0)
                Z1 = Z
                Z  = Z1 - P1 / PP
                if (abs(Z - Z1) < TOL) exit
            end do
            QSI(I)           = -Z
            QSI(NPG + 1 - I) =  Z
            W(I)             = 2.0 / ((1.0 - Z * Z) * PP * PP)
            W(NPG + 1 - I)   = W(I)
        end do
    end subroutine GAUSS_LEGENDRE

    !-------------------------------------------------------------------
    !> @brief Exponential integral function E1(AF) using Abramowitz &
    !>        Stegun polynomial/rational approximations in double precision.
    !-------------------------------------------------------------------
    pure real function EXPI(AF)
        implicit none
        real, intent(in) :: AF

        real, parameter :: C1  = 0.57721566,   C2  = 0.99999193
        real, parameter :: C3  = 0.24991055,   C4  = 0.05519968
        real, parameter :: C5  = 0.00976004,   C6  = 0.00107857
        real, parameter :: C7  = 8.5733287401, C8  = 18.0590169730
        real, parameter :: C9  = 8.6347608925, C10 = 0.2677737343
        real, parameter :: C11 = 9.5733223454, C12 = 25.6329561486
        real, parameter :: C13 = 21.0996530827, C14 = 3.9584969228

        if (AF > 199.0) then
            EXPI = 0.0
        else if (AF <= 1.0) then
            EXPI = -C1 + C2*AF - C3*AF**2 + C4*AF**3 - C5*AF**4 + C6*AF**5 - log(AF)
        else
            EXPI = (AF**4 + C7*AF**3 + C8*AF**2 + C9*AF + C10) / &
                   ((AF**4 + C11*AF**3 + C12*AF**2 + C13*AF + C14) * AF * exp(AF))
        end if
    end function EXPI

    !-------------------------------------------------------------------
    !> @brief Smooth part of E1(AF) + ln(AF) used for analytic singular
    !>        integration in the first time step.
    !-------------------------------------------------------------------
    pure real function E1S(AF)
        implicit none
        real, intent(in) :: AF

        real, parameter :: C1  = 0.57721566,   C2  = 0.99999193
        real, parameter :: C3  = 0.24991055,   C4  = 0.05519968
        real, parameter :: C5  = 0.00976004,   C6  = 0.00107857
        real, parameter :: C7  = 8.5733287401, C8  = 18.0590169730
        real, parameter :: C9  = 8.6347608925, C10 = 0.2677737343
        real, parameter :: C11 = 9.5733223454, C12 = 25.6329561486
        real, parameter :: C13 = 21.0996530827, C14 = 3.9584969228

        if (AF > 199.0) then
            E1S = 0.0
        else if (AF <= 1.0) then
            E1S = -C1 + C2*AF - C3*AF**2 + C4*AF**3 - C5*AF**4 + C6*AF**5
        else
            E1S = (AF**4 + C7*AF**3 + C8*AF**2 + C9*AF + C10) / &
                  ((AF**4 + C11*AF**3 + C12*AF**2 + C13*AF + C14) * AF * exp(AF)) &
                  + log(AF)
        end if
    end function E1S

    !-------------------------------------------------------------------
    !> @brief Hyperbolic sine (sinh) coordinate transformation
    !>        (Johnston & Elliott, 2005) for near-singular element integrals.
    !-------------------------------------------------------------------
    pure subroutine SINHTRANS(U, A_IN, B_IN, XI_OUT, JT)
        implicit none
        real, intent(in)  :: U, A_IN, B_IN
        real, intent(out) :: XI_OUT, JT
        real :: B, MI, ETA

        B = max(B_IN, 1.0e-8)
        MI     = 0.5 * (asinh((1.0 + A_IN) / B) + asinh((1.0 - A_IN) / B))
        ETA    = 0.5 * (asinh((1.0 + A_IN) / B) - asinh((1.0 - A_IN) / B))
        XI_OUT = A_IN + B * sinh(MI * U - ETA)
        JT     = B * MI * cosh(MI * U - ETA)
    end subroutine SINHTRANS

    !-------------------------------------------------------------------
    !> @brief Precomputes element midpoints, lengths, outward unit normals,
    !>        and discrete observation times TF(1:NINTER).
    !-------------------------------------------------------------------
    subroutine PRECOMPUTE_GEOMETRY()
        implicit none
        integer :: J, NO1, NO2, IT
        real :: DX, DY, L_J

        do J = 1, NE
            NO1 = CONEC(J, 1)
            NO2 = CONEC(J, 2)
            DX = XG(NO2) - XG(NO1)
            DY = YG(NO2) - YG(NO1)
            L_J = sqrt(DX*DX + DY*DY)
            XM(J)     = 0.5 * (XG(NO1) + XG(NO2))
            YM(J)     = 0.5 * (YG(NO1) + YG(NO2))
            EL_LEN(J) = L_J
            NX_VEC(J) =  DY / L_J
            NY_VEC(J) = -DX / L_J
        end do

        do IT = 1, NINTER
            TF(IT) = real(IT) * DELTA_T
        end do
    end subroutine PRECOMPUTE_GEOMETRY

    !-------------------------------------------------------------------
    !> @brief Computes transient BEM influence matrices G and H for all
    !>        time steps in parallel using OpenMP.
    !>        Geometry (midpoints, lengths, normals, Gauss points, and
    !>        source-field distances) is computed once per (J, I) pair and
    !>        reused across all NINTER time steps.
    !-------------------------------------------------------------------
    subroutine MATGH()
        implicit none
        real, allocatable :: QSI(:), W(:), XK_ALL(:,:), YK_ALL(:,:)
        real, allocatable :: RA(:), DN(:)
        real :: N1, N2, AE(2), TA(2), TAL
        real :: IU, IQ, G1, H1, G1R, RX, RY, HALF_LEN, R2
        integer :: I, J, IT, I1, KM, NO1, NO2
        logical :: SINGULAR

        allocate(QSI(NG), W(NG))
        allocate(XK_ALL(NG, NE), YK_ALL(NG, NE))

        call GAUSS_LEGENDRE(NG, QSI, W)

        ! Precompute Gauss integration point coordinates on all elements
        do J = 1, NE
            NO1 = CONEC(J, 1)
            NO2 = CONEC(J, 2)
            do I1 = 1, NG
                N1 = 0.5 * (1.0 - QSI(I1))
                N2 = 0.5 * (1.0 + QSI(I1))
                XK_ALL(I1, J) = N1 * XG(NO1) + N2 * XG(NO2)
                YK_ALL(I1, J) = N1 * YG(NO1) + N2 * YG(NO2)
            end do
        end do

        TA(1) = 0.0
        TA(2) = DELTA_T

        ! Parallelize assembly over (element J, collocation node I)
        !$omp parallel do collapse(2) default(shared) &
        !$omp private(J, I, IT, I1, KM, RX, RY, RA, DN, TAL, AE, &
        !$omp         SINGULAR, IU, IQ, G1, H1, G1R, HALF_LEN, R2) &
        !$omp schedule(dynamic)
        do J = 1, NE
            do I = 1, NE
                allocate(RA(NG), DN(NG))
                HALF_LEN = 0.5 * EL_LEN(J)

                ! Compute distance and normal derivative once for pair (I, J)
                do I1 = 1, NG
                    RX     = XK_ALL(I1, J) - XM(I)
                    RY     = YK_ALL(I1, J) - YM(I)
                    RA(I1) = sqrt(RX*RX + RY*RY)
                    DN(I1) = -(RX * NX_VEC(J) + RY * NY_VEC(J))
                end do

                ! Evaluate influence coefficients across all time steps
                do IT = 1, NINTER
                    SINGULAR = (I == J) .and. (IT == 1)
                    G1 = 0.0
                    H1 = 0.0

                    do I1 = 1, NG
                        do KM = 1, 2
                            TAL = TF(IT) - TA(KM)
                            if (K_DIFF * TAL < 1.0e-12) then
                                AE(KM) = 200.0
                            else
                                AE(KM) = (RA(I1)**2) / (4.0 * K_DIFF * TAL)
                            end if
                        end do

                        if (SINGULAR) then
                            IU = E1S(AE(1)) / (4.0 * PI)
                            IQ = 0.0
                        else
                            IU = (EXPI(AE(1)) - EXPI(AE(2))) / (4.0 * PI)
                            R2 = RA(I1)**2
                            if (R2 > 1.0e-24) then
                                IQ = DN(I1) * (exp(-AE(1)) - exp(-AE(2))) / (2.0 * PI * R2)
                            else
                                IQ = 0.0
                            end if
                        end if

                        G1 = G1 + IU * W(I1) * HALF_LEN
                        H1 = H1 + IQ * W(I1) * HALF_LEN
                    end do

                    if (SINGULAR) then
                        G1R = EL_LEN(J) * (log(HALF_LEN) - 1.0 &
                              - 0.5 * log(4.0 * K_DIFF * DELTA_T)) / (2.0 * PI)
                        G1  = G1 - G1R
                        H1  = 0.5
                    end if

                    G(I, J, IT) = G1
                    H(I, J, IT) = H1
                end do

                deallocate(RA, DN)
            end do
        end do
        !$omp end parallel do

        deallocate(QSI, W, XK_ALL, YK_ALL)
    end subroutine MATGH

    !-------------------------------------------------------------------
    !> @brief Forms the constant system matrix C_MAT (factored once via
    !>        LAPACK DGETRF) and marches through NINTER time steps solving
    !>        each step with LAPACK DGETRS.
    !-------------------------------------------------------------------
    subroutine SOLVE_TRANSIENT()
        implicit none
        integer :: I, J, IT, IT1, L

        ! Build C_MAT (unknown coefficients) and B_MAT (prescribed coefficients)
        do J = 1, NE
            if (KODE(J) == 0) then
                C_MAT(:, J) = -G(:, J, 1)
                B_MAT(:, J) = -H(:, J, 1)
            else
                C_MAT(:, J) =  H(:, J, 1)
                B_MAT(:, J) =  G(:, J, 1)
            end if
        end do

        ! Copy C_MAT into A_MAT and compute LU factorization ONCE
        A_MAT = C_MAT
        call FACTOR_SYSTEM(NE, A_MAT, NE, IPIV, IPR)

        ! Time-marching loop
        do IT = 1, NINTER
            do I = 1, NE
                if (KODE(I) == 0) then
                    VPR(I) = FI(I, IT)
                else
                    VPR(I) = DFI(I, IT)
                end if
            end do

            ! First term of RHS: B_MAT * VPR
            RHS = matmul(B_MAT, VPR)

            ! Accumulate history contributions from previous time steps
            if (IT > 1) then
                IT1 = IT + 1
                do L = 2, IT
                    RHS = RHS - matmul(H(:, :, L), FI(:, IT1 - L)) &
                              + matmul(G(:, :, L), DFI(:, IT1 - L))
                end do
            end if

            ! Solve factored linear system A * X = RHS using LAPACK DGETRS
            call SOLVE_FACTORED(NE, A_MAT, NE, IPIV, RHS, IPR)

            ! Store solved boundary unknowns
            do I = 1, NE
                if (KODE(I) == 0) then
                    DFI(I, IT) = RHS(I)
                else
                    FI(I, IT)  = RHS(I)
                end if
            end do
        end do
    end subroutine SOLVE_TRANSIENT

    !-------------------------------------------------------------------
    !> @brief Evaluates temperature UPI(1:NPI, 1:NINTER) at all internal
    !>        points in parallel using OpenMP, with sinh transformation
    !>        for points near the boundary.
    !-------------------------------------------------------------------
    subroutine PTINT()
        implicit none
        real, allocatable :: QSI(:), W(:)
        integer :: KN

        if (NPI <= 0) return

        allocate(QSI(NG), W(NG))
        call GAUSS_LEGENDRE(NG, QSI, W)

        !$omp parallel do default(shared) private(KN) schedule(dynamic)
        do KN = 1, NPI
            call EVAL_INTERNAL_POINT(KN, QSI, W)
        end do
        !$omp end parallel do

        deallocate(QSI, W)
    end subroutine PTINT

    !-------------------------------------------------------------------
    !> @brief Computes transient temperature history for a single internal
    !>        point KN across all time steps.
    !-------------------------------------------------------------------
    subroutine EVAL_INTERNAL_POINT(KN, QSI, W)
        implicit none
        integer, intent(in) :: KN
        real,    intent(in) :: QSI(NG), W(NG)

        real :: GPI_LOC(NE, NINTER), HPI_LOC(NE, NINTER)
        real :: TA(2), AE(2), TAL
        real :: XP, YP, X1, Y1, X2, Y2, DX, DY, L_J, L2, HALF_LEN
        real :: XI_PROJ, XI_CLAMP, PSX, PSY, XC, YC, D_PERP, D_SEG
        real :: A_SINH, B_SINH, XI_G, JT, W_EFF
        real :: X_INT, Y_INT, RX, RY, RA, R2, DN, IU, IQ, G1, H1
        integer :: J, IT, IT1, L, I1, KM, NO1, NO2
        logical :: USE_SINH

        XP = XI(KN)
        YP = YI(KN)
        TA(1) = 0.0
        TA(2) = DELTA_T

        UPI(KN, :) = 0.0

        do J = 1, NE
            NO1 = CONEC(J, 1)
            NO2 = CONEC(J, 2)
            X1 = XG(NO1); Y1 = YG(NO1)
            X2 = XG(NO2); Y2 = YG(NO2)
            DX = X2 - X1
            DY = Y2 - Y1
            L2 = DX*DX + DY*DY
            L_J = EL_LEN(J)
            HALF_LEN = 0.5 * L_J

            ! Distance from internal point to element J for near-singularity check
            XI_PROJ  = 2.0 * ((XP - X1) * DX + (YP - Y1) * DY) / L2 - 1.0
            XI_CLAMP = min(max(XI_PROJ, -1.0), 1.0)
            XC = 0.5 * (1.0 - XI_CLAMP) * X1 + 0.5 * (1.0 + XI_CLAMP) * X2
            YC = 0.5 * (1.0 - XI_CLAMP) * Y1 + 0.5 * (1.0 + XI_CLAMP) * Y2
            D_SEG = sqrt((XC - XP)**2 + (YC - YP)**2)

            USE_SINH = (D_SEG < 1.5 * L_J)
            if (USE_SINH) then
                if (abs(XI_PROJ) <= 1.2) then
                    PSX = 0.5 * (1.0 - XI_PROJ) * X1 + 0.5 * (1.0 + XI_PROJ) * X2
                    PSY = 0.5 * (1.0 - XI_PROJ) * Y1 + 0.5 * (1.0 + XI_PROJ) * Y2
                    D_PERP = sqrt((PSX - XP)**2 + (PSY - YP)**2)
                    A_SINH = XI_PROJ
                    B_SINH = max(2.0 * D_PERP / L_J, 1.0e-6)
                else
                    A_SINH = XI_CLAMP
                    B_SINH = max(2.0 * D_SEG / L_J, 1.0e-6)
                end if
            end if

            do IT = 1, NINTER
                G1 = 0.0
                H1 = 0.0

                do I1 = 1, NG
                    if (USE_SINH) then
                        call SINHTRANS(QSI(I1), A_SINH, B_SINH, XI_G, JT)
                        W_EFF = W(I1) * JT * HALF_LEN
                    else
                        XI_G  = QSI(I1)
                        W_EFF = W(I1) * HALF_LEN
                    end if

                    X_INT = 0.5 * (1.0 - XI_G) * X1 + 0.5 * (1.0 + XI_G) * X2
                    Y_INT = 0.5 * (1.0 - XI_G) * Y1 + 0.5 * (1.0 + XI_G) * Y2

                    RX = X_INT - XP
                    RY = Y_INT - YP
                    R2 = RX*RX + RY*RY
                    RA = sqrt(R2)
                    DN = -(RX * NX_VEC(J) + RY * NY_VEC(J))

                    do KM = 1, 2
                        TAL = TF(IT) - TA(KM)
                        if (K_DIFF * TAL < 1.0e-12) then
                            AE(KM) = 200.0
                        else
                            AE(KM) = R2 / (4.0 * K_DIFF * TAL)
                        end if
                    end do

                    IU = (EXPI(AE(1)) - EXPI(AE(2))) / (4.0 * PI)
                    if (R2 > 1.0e-24) then
                        IQ = DN * (exp(-AE(1)) - exp(-AE(2))) / (2.0 * PI * R2)
                    else
                        IQ = 0.0
                    end if

                    G1 = G1 + IU * W_EFF
                    H1 = H1 + IQ * W_EFF
                end do

                GPI_LOC(J, IT) = G1
                HPI_LOC(J, IT) = H1

                IT1 = IT + 1
                do L = 1, IT
                    UPI(KN, IT) = UPI(KN, IT) &
                                  - HPI_LOC(J, L) * FI(J, IT1 - L) &
                                  + GPI_LOC(J, L) * DFI(J, IT1 - L)
                end do
            end do
        end do
    end subroutine EVAL_INTERNAL_POINT

end module assembly_mod
