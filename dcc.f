C (* ***************************************************************** *)
C (*                                                                   *)
C (*   COPPE/UFRJ - PROGRAMA DE ENGENHARIA CIVIL                       *)
C (*                                                                   *)
C (*   Solução da Equação da Difusão pelo Método dos Elementos de Contorno,*)
C (*   usando elementos constantes no espaço e no tempo                *)
C (*                                                                   *)
C (*   Teresinha Costa Effren                                          *)
C (*                                                                   *)
C (* ***************************************************************** *)
C ======================================================================
      PROGRAM DCC
C ======================================================================
C
      INTEGER NDELTM          ! No. Máx. de Intervalos de Tempo
      INTEGER NELMX           ! No. Máx. de elementos
      INTEGER NGAMX           ! No. Máx. de pontos de Gauss
      INTEGER NNGM            ! No. Máx. de nós geométricos
      INTEGER NPINTM          ! No. Máx. de pontos internos
C
      REAL*8 DeltaT           ! Intervalo de tempo em segundos
C
      CHARACTER*10 FILEIN     ! Nome do arquivo de entrada
      CHARACTER*10 FILEOUT    ! Nome do arquivo de saída
C
      INTEGER Inic            ! Arquivo de inicialização
      INTEGER Inp             ! Arquivo de entrada
      INTEGER It              ! Intervalo de tempo atual
      INTEGER Ipr             ! Arquivo de saída
C
      REAL*8 K                ! Coeficiente
C
      INTEGER NDMX            ! No. máx. de intervalos de tempo
      INTEGER Ne              ! No. de elementos
      INTEGER NEMX            ! No. máx. de elementos
      INTEGER NG              ! No. de pontos de Gauss
      INTEGER NGMX            ! No. máx. de nós geométricos
      INTEGER Ninter          ! No. de intervalos tempo do arq.
      INTEGER Nn              ! No. de nós geométricos
      INTEGER Npi             ! No. de pontos internos
      INTEGER NPIMX           ! No. máx de pontos internos
C
      REAL*8 Pi               ! Valor de Pi
C
      PARAMETER (NDELTM = 40)
      PARAMETER (NELMX = 60)
      PARAMETER (NGAMX = 4)
      PARAMETER (NNGM = 61)
      PARAMETER (NPINTM = 10)
      PARAMETER (Pi = 3.141592654 D0)
C
      REAL*8 Ae(2)            ! Argumento da exponencial integral
      REAL*8 Cg(NGAMX)        ! Coordenadas de Gauss
      REAL*8 F(NELMX)         ! Vetor de termos indep. no tempo atual
      INTEGER Kode(NELMX)     ! Código do valor prescrito
      REAL*8 Len(NELMX)       ! Comprimento de elemento
      REAL*8 Ome(NGAMX)       ! Pesos de Gauss
      REAL*8 Tf(NDELTM)       ! Tempo de contagem da resposta
      REAL*8 Valor(NDELTM)    ! Leitura de valores prescritos
      REAL*8 Vpr(NELMX)       ! Val. prescritos no tempo atual
C
      REAL*8 Xg(NNGM), Yg(NNGM)     ! Coordenadas dos nós geométricos
      REAL*8 Xi(NPINTM), Yi(NPINTM) ! Coordenadas dos pontos internos
      REAL*8 Xk(NGAMX), Yk(NGAMX)   ! Coord. dos pontos de Gauss no elem.
      REAL*8 Xm(NELMX), Ym(NELMX)   ! Coord. do ponto medio do elem.
C
      REAL*8 A(NELMX,NELMX)         ! Matriz do sist. que entra no solver
      REAL*8 B(NELMX,NELMX)         ! Multiplica os valores prescritos
      REAL*8 C(NELMX,NELMX)         ! Back-up da matriz [A]
      REAL*8 Dfi(NELMX,NDELTM)      ! Valores da der. normal(P)
      REAL*8 Fi(NELMX,NDELTM)       ! Valores de potencial(U)
      REAL*8 Gpi(NELMX,NDELTM)      ! Coef.influência G no ponto int.
      REAL*8 Hpi(NELMX,NDELTM)      ! Coef.influência H no ponto int.
C
      REAL*8 G(NELMX,NELMX,NDELTM)  ! Matriz de coeficientes Gijt
      REAL*8 H(NELMX,NELMX,NDELTM)  ! Matriz de coeficientes Hijt
      INTEGER Inc(NELMX,2)          ! Matriz de conectividades
      REAL*8 Upi(NPINTM,NDELTM)     ! Potencial dos pontos internos
C
      NDMX = NDELTM
      NEMX = NELMX
      NPIMX = NPINTM
      NG = NGAMX
      NGMX = NNGM
C
C Leitura dos arquivos
C
      Inic= 4
      Inp = 5
      Ipr = 6
C
      OPEN(Inic, FILE = 'INICIO.DAT', STATUS = 'OLD')
      READ(Inic, '(A)') FILEIN
      READ(Inic, '(A)') FILEOUT
      OPEN(Inp, FILE = FILEIN, STATUS = 'OLD')
      CLOSE (Inic)
C
C LÊ OS DADOS DE ENTRADA
C
C ---------------------------------------------
      CALL LER (
     D     NDMX,NEMX,NGMX,NPIMX,
     S     DeltaT,Inp,Ipr,K,Ne,Ninter,Nn,Npi,
     V     Kode,
     V     Len,Valor,Xg,Xm,Xi,Yg,Ym,Yi,
     M     Dfi,Fi,Inc        )
C ---------------------------------------------
C
C VALORES DE COORDENADAS E PONTOS DE GAUSS
C
C ---------------------------------------------
      CALL GAUSS (
     D     NG,
     V     Cg, Ome)
C ----------------
C
      DO 5 It=1,NInter
       Tf(It)= It*DeltaT
C
C CÁLCULO DOS TERMOS DAS MATRIZES H E G EM CADA INTERVALO
C
C -----------------------------------------
      CALL MATGH (
     D     NDMX,NEMX,NG,NGMX,
     S     DeltaT,It,K,Ne,Pi,
     V     Ae,Cg,Len,Ome,
     V     Tf,Xg,Xk,Xm,Yg,Yk,Ym,
     M     G,H,Inc        )
C ------------------------------------------
C
    5 CONTINUE
C
C Preenchimento das matrizes C e B que multiplicam
C as incógnitas e os valores prescritos
C
      DO 20 j=1,Ne
       IF (Kode(j).EQ.0) THEN
        DO 10 i=1,Ne
         C(i,j)= -G(i,j,1)
         B(i,j)= -H(i,j,1)
   10   CONTINUE
       ELSE
        DO 15 i=1,Ne
         C(i,j)= H(i,j,1)
         B(i,j)= G(i,j,1)
   15   CONTINUE
       END IF
   20 CONTINUE

      DO 900 It=1,Ninter
C
C Preenchimento da matriz que multiplica as incógnitas
C e do vetor de condições de contorno prescritas(vpr)
C
       DO 200 I= 1,Ne
        DO 100 J= 1,Ne
         A(I,J)=C(I,J)
  100   CONTINUE
C
        IF (Kode(I).EQ.0)THEN
         Vpr(I)=Fi(I,It)
        ELSE
         Vpr(I)=Dfi(I,It)
        END IF
  200  CONTINUE
C
C Cálculo do 1o. termo do vetor de termos independentes
C B * VPR
C
       DO 400 I= 1,Ne
        F(i)= 0.0
        DO 300 J= 1,Ne
         F(i)=F(i) + B(i,j)*Vpr(j)
  300   CONTINUE
  400  CONTINUE
C
C Cálculo da contribuição do passado no vetor
C de termos independentes
C
       IF (It.GT.1)THEN
        DO 700 I= 1,Ne
         DO 600 J= 1,Ne
          It1=It + 1
          DO 500 L= 2,It
           F(I) = F(I) - H(I,J,L)*Fi(J,It1-L)
     *                 + G(I,J,L)*Dfi(J,It1-L)
  500     CONTINUE
  600    CONTINUE
  700   CONTINUE
       END IF
C
C SOLVER PARA RESOLUÇÃO DO SISTEMA DE EQUAÇÕES
C
C -------------------------
       CALL SLNPD (
     S      NEMX,Ne,
     V      A,F )
C -------------------------
C
C Preenchimento das matrizes de potencial(Fi) e fluxo(Dfi)
C
       DO 800 I=1,Ne
        IF (Kode(I).EQ.0)THEN
         Dfi(I,It)=F(I)
        ELSE
         Fi(I,It)=F(I)
        END IF
  800  CONTINUE
  900 CONTINUE
      IF (Npi.GT.0)THEN
C
C CALCULA O POTENCIAL DOS PONTOS INTERNOS
C
C --------------------------------------------------
       CALL PTINT (
     D      NDMX,NEMX,NG,NGMX,NPIMX,
     S      DeltaT,K,Ne,Ninter,Npi,Pi,
     V      Ae,Cg, Len, Ome,
     V      Tf,Xg,Xi,Xk,Yg,Yi,Yk,
     M      Dfi,Fi,Gpi,Hpi, Inc,Upi  )
C -----------------------------------------------
C
      END IF
C
      OPEN(Ipr, FILE = FILEOUT, STATUS = 'NEW')
C
C CRIAÇÃO DE UM ARQUIVO COM OS RESULTADOS
C
C ------------------------------------------
      CALL RESULT (
     D      NDMX,NEMX,NPIMX,
     S      Ipr,Ne,Ninter,Npi,
     V      Xi,Xm,Yi,Ym,
     M      Dfi,Fi,Upi    )
C -----------------------------------------
C
      CLOSE (Ipr)
      STOP
      END
C
C
C *********************************
      SUBROUTINE GAUSS (
     D      NG,
     V      Cg, Ome)
C *********************************
C
      INTEGER NG              ! Número de pontos de Gauss
      REAL*8 Cg(NG)           ! Coordenadas de Gauss
      REAL*8 Ome(NG)          ! Pesos de Gauss
C
      Cg(1)= 0.86113631159 D0
      Cg(2)= -Cg(1)
      Cg(3)= 0.33998104358 D0
      Cg(4)= -Cg(3)
      Ome(1)= 0.34785484514 D0
      Ome(2)= Ome(1)
      Ome(3)= 0.65214515486 D0
      Ome(4)= Ome(3)
C
      RETURN
      END
C
C
C *******************************************************
      FUNCTION EXPI (
     V      Af )
C *******************************************************
C
C CALCULA A FUNÇÃO EXPONENCIAL INTEGRAL
C
      REAL*8 Af
      REAL*8 C1,C2,C3,C4,C5,C6,C7,C8
      REAL*8 C9,C10,C11,C12,C13,C14
C
      PARAMETER (C1 = 0.57721566D0, C2 = 0.99999193D0)
      PARAMETER (C3 = 0.24991055D0, C4 = 0.05519968D0)
      PARAMETER (C5 = 0.00976004D0, C6 = 0.00107857D0)
C
      PARAMETER (C7 = 8.5733287401D0, C8 = 18.0590169730D0)
      PARAMETER (C9 = 8.6347608925D0, C10 = 0.2677737343D0)
      PARAMETER (C11 = 9.5733223454D0, C12 = 25.6329561486D0)
      PARAMETER (C13 = 21.0996530827D0, C14 = 3.9584969228D0)
C
      IF (Af.GT.199) THEN
       ExpI=0.0D0
      ELSE
       IF(Af.LE.1) THEN
        ExpI=-C1+C2*Af-C3*Af**2+C4*Af**3-C5*Af**4+C6*Af**5-DLOG(Af)
       ELSE
        ExpI=(Af**4+C7*Af**3+C8*Af**2+C9*Af+C10)/
     /       ((Af**4+C11*Af**3+C12*Af**2+C13*Af+C14)*Af*EXP(Af))
       END IF
      END IF
C
      RETURN
      END
C
C *******************************************************
      FUNCTION E1S(
     S      Af )
C *******************************************************
C
C CALCULA A FUNÇÃO EXPONENCIAL INTEGRAL SOMADA AO
C LN(INTERVALO 1)
C
      REAL*8 Af
      REAL*8 C1,C2,C3,C4,C5,C6,C7,C8
      REAL*8 C9,C10,C11,C12,C13,C14
C
      PARAMETER (C1 = 0.57721566D0, C2 = 0.99999193D0)
      PARAMETER (C3 = 0.24991055D0, C4 = 0.05519968D0)
      PARAMETER (C5 = 0.00976004D0, C6 = 0.00107857D0)
C
      PARAMETER (C7 = 8.5733287401D0, C8 = 18.0590169730D0)
      PARAMETER (C9 = 8.6347608925D0, C10 = 0.2677737343D0)
      PARAMETER (C11 = 9.5733223454D0, C12 = 25.6329561486D0)
      PARAMETER (C13 = 21.0996530827D0, C14 = 3.9584969228D0)
C
      IF (Af.GT.199) THEN
       E1S=0.0D0
      ELSE
       IF(Af.LE.1) THEN
        E1S=-C1+C2*Af-C3*Af**2+C4*Af**3-C5*Af**4+C6*Af**5
       ELSE
        E1S=(Af**4+C7*Af**3+C8*Af**2+C9*Af+C10)/
     /       ((Af**4+C11*Af**3+C12*Af**2+C13*Af+C14)*Af*EXP(Af))
     +       + DLOG(Af)
       END IF
      END IF
C
      RETURN
      END
C
C *****************************************************
      SUBROUTINE LER (
     D      NDMX,NEMX,NGMX,NPIMX,
     S      DeltaT,Inp,Ipr,K,Ne,Ninter,Nn,Npi,
     V      Kode,
     V      Len,Valor,Xg,Xm,Xi,Yg,Ym,Yi,
     M      Dfi,Fi,Inc        )
C *****************************************************
C
C ROTINA PARA LEITURA DOS DADOS DE ENTRADA
C
      INTEGER NDMX, NEMX, NPIMX, NGMX
      CHARACTER*80 TITLE
      REAL*8 DeltaT,K
      INTEGER Ke, Kn, Inp,Ipr,It
      INTEGER Ne,Ninter,Nn,No1,No2,Npi
      REAL*8 Xg1,Xg2,Yg1,Yg2
C
      INTEGER Kode(NEMX)
      REAL*8 Len(NEMX)
      REAL*8 Valor(NDMX)
      REAL*8 Xg(NGMX),Yg(NGMX)
      REAL*8 Xi(NPIMX),Yi(NPIMX)
      REAL*8 Xm(NEMX),Ym(NEMX)
      REAL*8 Dfi(NEMX,NDMX),Fi(NEMX,NDMX)
      INTEGER Inc(NEMX,2)
C
C Lê e imprime título
C
      READ(Inp,'(A)') TITLE
      WRITE(Ipr,'(A)') TITLE
C
      READ(Inp,*) K
      READ(Inp,*) DeltaT
      READ(Inp,*) Nn
C
      DO 100 I = 1,Nn
       READ(Inp,*) Kn, Xg(Kn), Yg(Kn)
  100 CONTINUE
C
      READ(Inp,*) Ne
      DO 200 I= 1,Ne
       READ(Inp,*) Ke, Inc(Ke,1),Inc(Ke,2),Kode(Ke)
       No1 = Inc(Ke,1)
       No2 = Inc(Ke,2)
       Xg1 = Xg(No1)
       Xg2 = Xg(No2)
       Yg1 = Yg(No1)
       Yg2 = Yg(no2)
       Xm(Ke) = (Xg1 + Xg2) / 2.D0
       Ym(Ke) = (Yg1 + Yg2) / 2.D0
       Len(Ke) = SQRT((Xg2-Xg1)**2 + (Yg2-Yg1)**2)
  200 CONTINUE
C
      READ(Inp,*) Ninter
      DO 500 I = 1,Ne
       READ(Inp,*) (Valor(It), It=1,NInter)
       IF (Kode(I).EQ.0)THEN
        DO 300 It = 1, Ninter
         Fi(I,It)= Valor(It)
  300   CONTINUE
       ELSE
        DO 400 It = 1,Ninter
         Dfi(I,It)= Valor(It)
  400   CONTINUE
       END IF
  500 CONTINUE
C
      READ(INP,*) Npi
      DO 600 I = 1,Npi
       READ(INP,*) Kn, Xi(Kn), Yi(Kn)
  600 CONTINUE
C
      CLOSE(Inp)
      RETURN
      END
C
C *****************************************************
      SUBROUTINE MATGH (
     D      NDMX,NEMX,NG,NGMX,
     S      DeltaT,It,K,Ne,Pi,
     V      Ae,Cg,Len,Ome,
     V      Tf,Xg,Xk,Xm,Yg,Yk,Ym,
     M      G,H, Inc     )
C *****************************************************
C
C ROTINA PARA CÁLCULO DOS ELEMENTOS DAS MATRIZES H E G
C
      INTEGER NDMX, NEMX, NG, NGMX
      INTEGER I,I1,It,J,Km,Ne, No1, No2
      REAL*8 Dn,DeltaT,G1R,G1,H1,IU,IQ, K,Pi, Ra,Tal
      LOGICAL Singular
C
      REAL*8 Ae(2)            ! Argumento de entrada em ExpI e E1S
      REAL*8 Cg(NG),Ome(NG)   ! Valor e peso dos pontos de Gauss
      REAL*8 Len(NEMX)        ! Comprimento do elemento
      REAL*8 Ta(2)            ! Tempo de injeção
      REAL*8 Tf(NDMX)         ! Tempo de contagem da resposta
      REAL*8 Xg(NGMX),Yg(NGMX) ! Coord. dos nos geométricos
      REAL*8 Xk(NG), Yk(NG)   ! Coord. dos pontos de Gauss no elem.
      REAL*8 Xm(NEMX),Ym(NEMX) ! Coord. do ponto medio do elem.
C
      REAL*8 G(NEMX,NEMX,NDMX),H(NEMX,NEMX,NDMX)
      INTEGER Inc(NEMX,2)
C
C Inicio do calculo das matrizes
C
      DO 200 I= 1,Ne
       DO 100 J= 1,Ne
        Singular= (I.EQ.J).AND.(It.EQ.1)
C
C INTEGRAÇÃO PELA QUADRATURA DE GAUSS
C
        No1= Inc(J,1)
        No2= Inc(J,2)
        G1 = 0.0 D0
        H1 = 0.0 D0
C
        DO 50 I1=1,NG
         Xk(I1)= Xg(No1)*(1-Cg(I1))/2.0D0 +
     +           Xg(No2)*(1+Cg(I1))/2.0D0
         Yk(I1)= Yg(No1)*(1-Cg(I1))/2.0D0 +
     +           Yg(No2)*(1+Cg(I1))/2.0D0
         Ra= SQRT((Xk(I1)-Xm(I))**2+(Yk(I1)-Ym(I))**2)
         Dn= ((Xm(I)-Xk(I1))*(Yg(No2)-Yg(No1)) -
     -        (Ym(I)-Yk(I1))*(Xg(No2)-Xg(No1)))/Len(J)
C
C Calcula as integrais de u* e q* no tempo
C
         Ta(1)= 0.0 D0
         Ta(2)= DeltaT
         DO 25 Km = 1,2
          Tal= Tf(It) - Ta(Km)
C
C Caso Tal < 0.0 , ainda não houve efeito no ponto em análise
C
          IF ((K*Tal).LT.(0.000001)) THEN
           Ae(Km)= 200
          ELSE
           Ae(Km)= Ra**2/(4.0D0*K*Tal)
          END IF
   25    CONTINUE
C
C Calcula as integrais da solução e do fluxo fundamental
C
         IF (Singular) THEN
          IU= E1S(Ae(1))/(4.0D0*K*Pi)
          IQ= 0.0D0
         ELSE
          IU= (EXPI(Ae(1))-EXPI(Ae(2)))/(4.0D0*K*Pi)
          IQ= Dn * (EXP(-Ae(1))-EXP(-Ae(2)))/(2.0D0*Pi*K*Ra**2)
         END IF
C
         G1= G1 + IU *Ome(I1)*Len(J)/2.0D0
         H1= H1 + IQ *Ome(I1)*Len(J)/2.0D0
   50   CONTINUE
C
        IF (Singular) THEN
         G1r=Len(J)*(DLOG(Len(J)/2.0D0)-1-DLOG(4.D0*K*DeltaT)/2.D0)
     /        /(2*Pi*K)
         G1=G1-G1r
         H1=0.5D0
        END IF
        G(i,j,it)= G1
        H(i,j,it)= H1
  100  CONTINUE
  200 CONTINUE
C
      RETURN
      END
C
C *****************************************************
      SUBROUTINE PTINT (
     D      NDMX,NEMX,NG,NGMX,NPIMX,
     S      DeltaT,K,Ne,Ninter,Npi,Pi,
     V      Ae,Cg, Len, Ome,
     V      Tf,Xg,Xi,Xk,Yg,Yi,Yk,
     M      Dfi,Fi,Gpi,Hpi, Inc,Upi  )
C *****************************************************
C
C ROTINA PARA CÁLCULO DO POTENCIAL DOS PONTOS INTERNOS
C
      INTEGER NDMX, NEMX, NG, NGMX, NPIMX
      INTEGER I1,It,It1,J,Km,Kn,L
      INTEGER Ne,Ninter,No1,No2,Npi
      REAL*8 Dn,DeltaT,G1,H1,IU,IQ
      REAL*8 K,Pi, Ra,Tal
C
      REAL*8 Ae(2)            ! Argumento de entrada em ExpI e E1S
      REAL*8 Cg(NG),Ome(NG)   ! Valor e peso dos pontos de Gauss
      REAL*8 Len(NEMX)        ! Comprimento do elemento
      REAL*8 Ta(2)            ! Tempo de injeção
      REAL*8 Tf(NDMX)         ! Tempo de contagem da resposta
      REAL*8 Xg(NGMX),Yg(NGMX) ! Coord. dos nos geométricos
      REAL*8 Xi(NPIMX),Yi(NPIMX) ! Coord. dos pontos internos
      REAL*8 Xk(NG), Yk(NG)   ! Coord. dos pontos de Gauss no elem.
C
      REAL*8 Dfi(NEMX,NDMX),Fi(NEMX,NDMX)
      REAL*8 Gpi(NEMX,NDMX), Hpi(NEMX,NDMX)
      INTEGER Inc(NEMX,2)
      REAL*8 Upi(NPIMX,NDMX)  ! Potencial dos pontos internos
C
C Início do cálculo do potencial dos pontos internos
C
      DO 400 Kn=1,Npi
       DO 300 It=1,Ninter
        Upi(Kn,It)= 0.0D0
        DO 200 J= 1,Ne
C
C INTEGRAÇÃO PELA QUADRATURA DE GAUSS
C
         No1= Inc(J,1)
         No2= Inc(J,2)
         G1 = 0.0 D0
         H1 = 0.0 D0
C
         DO 50 I1=1,NG
          Xk(I1)= Xg(No1)*(1-Cg(I1))/2.0D0 +
     +            Xg(No2)*(1+Cg(I1))/2.0D0
          Yk(I1)= Yg(No1)*(1-Cg(I1))/2.0D0 +
     +            Yg(No2)*(1+Cg(I1))/2.0D0
          Ra= SQRT((Xk(I1)-Xi(Kn))**2+(Yk(I1)-Yi(Kn))**2)
          Dn= ((Xi(Kn)-Xk(I1))*(Yg(No2)-Yg(No1)) -
     -         (Yi(Kn)-Yk(I1))*(Xg(No2)-Xg(No1)))/Len(J)
C
C Calcula as integrais de u* e q* no tempo
C
          Ta(1)= 0.0 D0
          Ta(2)= DeltaT
          DO 25 Km = 1,2
           Tal= Tf(It) - Ta(Km)
C
C Caso Tal < 0.0 , ainda não houve efeito no ponto em análise
C
           IF ((K*Tal).LT.(0.000001)) THEN
            Ae(Km)= 200
           ELSE
            Ae(Km)= Ra**2/(4*K*Tal)
           END IF
   25     CONTINUE
C
C Calcula as integrais da solução e do fluxo fundamental
C
          IU= (EXPI(Ae(1))-EXPI(Ae(2)))/(4.0D0*K*Pi)
          IQ= Dn * (EXP(-Ae(1))-EXP(-Ae(2)))/(2.0D0*Pi*K*Ra**2)
C
          G1= G1 + IU *Ome(I1)*Len(J)/2.0D0
          H1= H1 + IQ *Ome(I1)*Len(J)/2.0D0
   50    CONTINUE
C
         Gpi(j,it)= G1
         Hpi(j,it)= H1
         It1= It + 1
         DO 100 L= 1,It
          Upi(Kn,It)= Upi(Kn,It) - Hpi(J,L)*Fi(J,It1-L)
     +                           + Gpi(J,L)*Dfi(J,It1-L)
  100    CONTINUE
  200   CONTINUE
  300  CONTINUE
  400 CONTINUE
C
      RETURN
      END
C
C SUBROTINA SLNPD. ESTA SUBROTINA É UMA ROTINA PADRÃO PARA
C RESOLVER MATRIZES NÃO POSITIVAS DEFINIDAS USANDO
C ELIMINAÇÃO DE GAUSS . SE A MATRIZ TIVER UM ZERO NA DIAGONAL
C A ROTINA TROCARÁ AS LINHAS, E SOMENTE RESOLVERÁ QUE O
C SISTEMA DA MATRIZ É SINGULAR QUANDO NENHUMA TROCA DE
C LINHA PRODUZA UM COEFICIENTE DIFERENTE DE ZERO.
C O RESULTADO SERÁ UM VETOR COM OS POTENCIAIS E DERIVADAS DO
C CONTORNO E SERÁ SOBRE ESCRITO EM DFI.
C
      SUBROUTINE SLNPD (
     S      MX,NN,
     V      A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C A = MATRIZ DO SISTEMA
C B = ORIGINALMENTE CONTÉM OS COEFICIENTES INDEPENDENTES.
C DEPOIS DA SOLUÇÃO CONTERÁ OS VALORES DESCONHECIDOS DO
C SISTEMA.
C NN = NÚMERO REAL DE INCÓGNITAS
C MX = DIMENSÃO DAS COLUNAS E LINHAS DE A
C
      DIMENSION A(MX,MX),B(MX)
C
      N1 = NN - 1
      DO 100 K = 1,N1
      K1 = K+1
      C = A(K,K)
      IF(ABS(C) - 0.000001)1,1,3
    1 DO 7 J = K1,NN
C
C MUDANÇAS DE LINHAS ATÉ ACHAR UM COEFICIENTE NA DIAGONAL
C DIFERENTE DE ZERO
C
      IF(ABS(A(J,K))-0.000001)7,7,5
    5 DO 6 L = K,NN
      C = A(K,L)
      A(K,L) = A(J,L)
    6 A(J,L) = C
      C = B(K)
      B(K) = B(J)
      B(J) = C
      C = A(K,K)
      GO TO 3
    7 CONTINUE
      GO TO 8
C
C DIVISÃO DA LINHA PELO COEFICIENTE DA DIAGONAL
C
    3 C = A(K,K)
      DO 4 J = K1,NN
    4 A(K,J) = A(K,J)/C
      B(K) = B(K)/C
C
C ELIMINAÇÃO DA INCÓGNITA DA LINHA I
C
      DO 10 I = K1,NN
      C = A(I,K)
      DO 9 J = K1,NN
    9 A(I,J) = A(I,J) - C*A(K,J)
   10 B(I) = B(I) - C*B(K)
  100 CONTINUE
C
C CÁLCULO DA ÚLTIMA INCÓGNITA
C
      IF(ABS(A(NN,NN))- 0.000001)8,8,101
  101 B(NN) = B(NN)/A(NN,NN)
C
C CÁLCULO DAS DEMAIS INCÓGNITAS PELO PROCESSO DE
C RETROSUBSTITUIÇÃO
C
      DO 200 L = 1,N1
      K = NN - L
      K1 = K + 1
      DO 200 J = K1,NN
  200 B(K) = B(K) - A(K,J)*B(J)
      GO TO 300
    8 WRITE (*,2) K
      WRITE (5,2) K
    2 FORMAT ('***** SINGULARIDADE NA COLUNA',I5)
  300 RETURN
      END
C
C *****************************************************
      SUBROUTINE RESULT (
     D      NDMX,NEMX,NPIMX,
     S      Ipr,Ne,Ninter,Npi,
     V      Xi,Xm,Yi,Ym,
     M      Dfi,Fi,Upi    )
C *****************************************************
C
C ROTINA PARA CRIAÇÃO DE UM ARQUIVO DE RESULTADOS
C
      INTEGER NDMX,NEMX,NPIMX
      INTEGER I,Ipr,It,Ne,Ninter,Npi
C
      REAL*8 Xi(NPIMX),Yi(NPIMX)
      REAL*8 Xm(NEMX),Ym(NEMX)
      REAL*8 Dfi(NEMX,NDMX),Fi(NEMX,NDMX)
      REAL*8 Upi(NPIMX,NDMX)
C
      DO 300 It= 1,Ninter
       WRITE(Ipr,*)'Tempo: ',It
       WRITE(Ipr,*) ' Elem.   Xm       Ym       U        P'
       DO 100 I= 1,Ne
        WRITE(Ipr,50) I,Xm(I),Ym(I),Fi(I,It),Dfi(I,It)
   50   FORMAT(1x,I4,3(2x,f10.3),2x,f11.4)
  100  CONTINUE
       IF (Npi.GT.0) THEN
        WRITE(Ipr,*) 'PInter.  Xi       Yi       U'
        DO 200 Kn= 1,Npi
         WRITE(Ipr,150) Kn,Xi(Kn),Yi(Kn),Upi(Kn,It)
  150    FORMAT(1x,I4,3(2x,f10.3))
  200   CONTINUE
       END IF
  300 CONTINUE
C
      RETURN
      END
