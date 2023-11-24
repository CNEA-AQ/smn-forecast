module megvea
   implicit none

   !Number of emission classes
   INTEGER, PARAMETER :: NCLASS = 19
   INTEGER, PARAMETER :: NEMIS  = NCLASS
   ! number of emission classes

   ! CO2 related emission activity factor parameters
   REAL,PARAMETER :: CO2   = 400.0
   REAL,PARAMETER :: ISmax =   1.344
   REAL,PARAMETER :: CO2h  =   1.4614
   REAL,PARAMETER :: Cstar = 585.0

   ! PSTD
   REAL,PARAMETER :: PSTD = 200

   ! canopy depth emission response
   REAL,PARAMETER :: CCD1 = -0.2
   REAL,PARAMETER :: CCD2 =  1.3

   !Light and temperature emission activity response coefficients for each emission class
   !LDF: light dependent fraction
   REAL           LDF(NCLASS)
   !CT1: temperature coefficient (emission type 1: light dependent)
   REAL           CT1(NCLASS)
   !Cleo: temperature coefficient (emission type 1: light dependent)
   REAL           Cleo(NCLASS)
   !beta: temperature coefficient (emission type 2: light independent)
   REAL           beta(NCLASS)

      DATA    beta(1),LDF(1),CT1(1),Cleo(1)        / 0.13,1.0,95,2    /
      DATA    beta(2),LDF(2),CT1(2),Cleo(2)        / 0.13,1.0,95,2    /
      DATA    beta(3),LDF(3),CT1(3),Cleo(3)        / 0.10,0.6,80,1.83 /
      DATA    beta(4),LDF(4),CT1(4),Cleo(4)        / 0.10,0.9,80,1.83 /
      DATA    beta(5),LDF(5),CT1(5),Cleo(5)        / 0.10,0.2,80,1.83 /
      DATA    beta(6),LDF(6),CT1(6),Cleo(6)        / 0.10,0.4,80,1.83 /
      DATA    beta(7),LDF(7),CT1(7),Cleo(7)        / 0.10,0.1,80,1.83 /
      DATA    beta(8),LDF(8),CT1(8),Cleo(8)        / 0.10,0.0,80,1.83 /
      DATA    beta(9),LDF(9),CT1(9),Cleo(9)        / 0.17,0.5,130,2.37/
      DATA    beta(10),LDF(10),CT1(10),Cleo(10)    / 0.17,0.4,130,2.37/
      DATA    beta(11),LDF(11),CT1(11),Cleo(11)    / 0.08,0.8,60,1.6  /
      DATA    beta(12),LDF(12),CT1(12),Cleo(12)    / 0.10,0.2,80,1.83 /
      DATA    beta(13),LDF(13),CT1(13),Cleo(13)    / 0.13,0.8,95,2    /
      DATA    beta(14),LDF(14),CT1(14),Cleo(14)    / 0.13,0.8,95,2    /
      DATA    beta(15),LDF(15),CT1(15),Cleo(15)    / 0.10,0.2,80,1.83 /
      DATA    beta(16),LDF(16),CT1(16),Cleo(16)    / 0.10,0.2,80,1.83 /
      DATA    beta(17),LDF(17),CT1(17),Cleo(17)    / 0.10,0.8,80,1.83 /
      DATA    beta(18),LDF(18),CT1(18),Cleo(18)    / 0.10,0.1,80,1.83 /
      DATA    beta(19),LDF(19),CT1(19),Cleo(19)    / 0.08,1.0,60,1.6  /

   ! Parameters for leaf age algorithm for each emission activity classes
   REAL           Anew(NCLASS)
   REAL           Agro(NCLASS)
   REAL           Amat(NCLASS)
   REAL           Aold(NCLASS)
      DATA Anew( 1), Agro( 1), Amat( 1), Aold( 1) / 0.05 ,0.6, 1.0, 0.9  /
      DATA Anew( 2), Agro( 2), Amat( 2), Aold( 2) / 0.05 ,0.6, 1.0, 0.9  /
      DATA Anew( 3), Agro( 3), Amat( 3), Aold( 3) / 2.0  ,1.8, 1.0, 1.05 /
      DATA Anew( 4), Agro( 4), Amat( 4), Aold( 4) / 2.0  ,1.8, 1.0, 1.05 /
      DATA Anew( 5), Agro( 5), Amat( 5), Aold( 5) / 2.0  ,1.8, 1.0, 1.05 /
      DATA Anew( 6), Agro( 6), Amat( 6), Aold( 6) / 2.0  ,1.8, 1.0, 1.05 /
      DATA Anew( 7), Agro( 7), Amat( 7), Aold( 7) / 2.0  ,1.8, 1.0, 1.05 /
      DATA Anew( 8), Agro( 8), Amat( 8), Aold( 8) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew( 9), Agro( 9), Amat( 9), Aold( 9) / 0.4  ,0.6, 1.0, 0.95 /
      DATA Anew(10), Agro(10), Amat(10), Aold(10) / 0.4  ,0.6, 1.0, 0.95 /
      DATA Anew(11), Agro(11), Amat(11), Aold(11) / 3.5  ,3.0, 1.0, 1.2  /
      DATA Anew(12), Agro(12), Amat(12), Aold(12) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(13), Agro(13), Amat(13), Aold(13) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(14), Agro(14), Amat(14), Aold(14) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(15), Agro(15), Amat(15), Aold(15) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(16), Agro(16), Amat(16), Aold(16) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(17), Agro(17), Amat(17), Aold(17) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(18), Agro(18), Amat(18), Aold(18) / 1.0  ,1.0, 1.0, 1.0  /
      DATA Anew(19), Agro(19), Amat(19), Aold(19) / 1.0  ,1.0, 1.0, 1.0  /

   !stress emission activity response coefficients for each emission class
   !CAQ: coefficient for poor Air Quality stress
   REAL           CAQ(NCLASS)
   !CHW: coefficient for high wind speed stress
   REAL           CHW(NCLASS)
   !CHT: coefficient for high temperature stress
   REAL           CHT(NCLASS)
   !CLT: coefficient for high temperature stress
   REAL           CLT(NCLASS)
      DATA    CAQ(1),CHW(1),CHT(1),CLT(1)           /  1,1,1,1  /
      DATA    CAQ(2),CHW(2),CHT(2),CLT(2)           /  1,1,1,1  /
      DATA    CAQ(3),CHW(3),CHT(3),CLT(3)           /  1,5,1,1  /
      DATA    CAQ(4),CHW(4),CHT(4),CLT(4)           /  5,5,5,5  /
      DATA    CAQ(5),CHW(5),CHT(5),CLT(5)           /  1,5,1,1  /
      DATA    CAQ(6),CHW(6),CHT(6),CLT(6)           /  1,5,1,1  /
      DATA    CAQ(7),CHW(7),CHT(7),CLT(7)           /  1,5,1,1  /
      DATA    CAQ(8),CHW(8),CHT(8),CLT(8)           /  1,1,1,1  /
      DATA    CAQ(9),CHW(9),CHT(9),CLT(9)           /  5,5,5,5  /
      DATA    CAQ(10),CHW(10),CHT(10),CLT(10)       /  5,5,5,5  /
      DATA    CAQ(11),CHW(11),CHT(11),CLT(11)       /  1,1,1,1  /
      DATA    CAQ(12),CHW(12),CHT(12),CLT(12)       /  1,1,1,1  /
      DATA    CAQ(13),CHW(13),CHT(13),CLT(13)       /  1,1,1,1  /
      DATA    CAQ(14),CHW(14),CHT(14),CLT(14)       /  1,1,1,1  /
      DATA    CAQ(15),CHW(15),CHT(15),CLT(15)       /  1,1,1,1  /
      DATA    CAQ(16),CHW(16),CHT(16),CLT(16)       /  1,1,1,1  /
      DATA    CAQ(17),CHW(17),CHT(17),CLT(17)       /  5,5,5,5  /
      DATA    CAQ(18),CHW(18),CHT(18),CLT(18)       /  1,1,1,1  /
      DATA    CAQ(19),CHW(19),CHT(19),CLT(19)       /  1,1,1,1  /

   !TAQ: threshold for poor Air Quality stress (ppm-hours)
   REAL           TAQ(NCLASS)
   !THW: threshold for high wind speed stress (m/s)
   REAL           THW(NCLASS)
   !THT: threshold for high temperature stress (Celsius degree)
   REAL           THT(NCLASS)
   !TLT: threshold for high temperature stress (Celsius degree)
   REAL           TLT(NCLASS)
      DATA    TAQ(1),THW(1),THT(1),TLT(1)           /  20,12,40,10  /
      DATA    TAQ(2),THW(2),THT(2),TLT(2)           /  20,12,40,10  /
      DATA    TAQ(3),THW(3),THT(3),TLT(3)           /  20,12,40,10  /
      DATA    TAQ(4),THW(4),THT(4),TLT(4)           /  20,12,40,10  /
      DATA    TAQ(5),THW(5),THT(5),TLT(5)           /  20,12,40,10  /
      DATA    TAQ(6),THW(6),THT(6),TLT(6)           /  20,12,40,10  /
      DATA    TAQ(7),THW(7),THT(7),TLT(7)           /  20,12,40,10  /
      DATA    TAQ(8),THW(8),THT(8),TLT(8)           /  20,12,40,10  /
      DATA    TAQ(9),THW(9),THT(9),TLT(9)           /  20,12,40,10  /
      DATA    TAQ(10),THW(10),THT(10),TLT(10)       /  20,12,40,10  /
      DATA    TAQ(11),THW(11),THT(11),TLT(11)       /  20,12,40,10  /
      DATA    TAQ(12),THW(12),THT(12),TLT(12)       /  20,12,40,10  /
      DATA    TAQ(13),THW(13),THT(13),TLT(13)       /  20,12,40,10  /
      DATA    TAQ(14),THW(14),THT(14),TLT(14)       /  20,12,40,10  /
      DATA    TAQ(15),THW(15),THT(15),TLT(15)       /  20,12,40,10  /
      DATA    TAQ(16),THW(16),THT(16),TLT(16)       /  20,12,40,10  /
      DATA    TAQ(17),THW(17),THT(17),TLT(17)       /  20,12,40,10  /
      DATA    TAQ(18),THW(18),THT(18),TLT(18)       /  20,12,40,10  /
      DATA    TAQ(19),THW(19),THT(19),TLT(19)       /  20,12,40,10  /

    !stress emission activity delta thresholds for each emission class
    !DTAQ: delta threshold for poor Air Quality stress (ppm-hours)
    REAL           DTAQ(NCLASS)
    !DTHW: delta threshold for high wind speed stress (m/s)
    REAL           DTHW(NCLASS)
    !DTHT: delta threshold for high temperature stress (Celsius degree)
    REAL           DTHT(NCLASS)
    !DTLT: delta threshold for low temperature stress (Celsius degree)
    REAL           DTLT(NCLASS)
      DATA    DTAQ(1),DTHW(1),DTHT(1),DTLT(1)               /  30,8,8,8  /
      DATA    DTAQ(2),DTHW(2),DTHT(2),DTLT(2)               /  30,8,8,8  /
      DATA    DTAQ(3),DTHW(3),DTHT(3),DTLT(3)               /  30,8,8,8  /
      DATA    DTAQ(4),DTHW(4),DTHT(4),DTLT(4)               /  30,8,8,8  /
      DATA    DTAQ(5),DTHW(5),DTHT(5),DTLT(5)               /  30,8,8,8  /
      DATA    DTAQ(6),DTHW(6),DTHT(6),DTLT(6)               /  30,8,8,8  /
      DATA    DTAQ(7),DTHW(7),DTHT(7),DTLT(7)               /  30,8,8,8  /
      DATA    DTAQ(8),DTHW(8),DTHT(8),DTLT(8)               /  30,8,8,8  /
      DATA    DTAQ(9),DTHW(9),DTHT(9),DTLT(9)               /  30,8,8,8  /
      DATA    DTAQ(10),DTHW(10),DTHT(10),DTLT(10)           /  30,8,8,8  /
      DATA    DTAQ(11),DTHW(11),DTHT(11),DTLT(11)           /  30,8,8,8  /
      DATA    DTAQ(12),DTHW(12),DTHT(12),DTLT(12)           /  30,8,8,8  /
      DATA    DTAQ(13),DTHW(13),DTHT(13),DTLT(13)           /  30,8,8,8  /
      DATA    DTAQ(14),DTHW(14),DTHT(14),DTLT(14)           /  30,8,8,8  /
      DATA    DTAQ(15),DTHW(15),DTHT(15),DTLT(15)           /  30,8,8,8  /
      DATA    DTAQ(16),DTHW(16),DTHT(16),DTLT(16)           /  30,8,8,8  /
      DATA    DTAQ(17),DTHW(17),DTHT(17),DTLT(17)           /  30,8,8,8  /
      DATA    DTAQ(18),DTHW(18),DTHT(18),DTLT(18)           /  30,8,8,8  /
      DATA    DTAQ(19),DTHW(19),DTHT(19),DTLT(19)           /  30,8,8,8  /





















contains

SUBROUTINE MEGVEA(  LAYERS, JDATE, ZTIME,                &
                    LAIp, LAIc,LDF_in,                              &
                                        GAMSM_in, MaxT, MinT,       &
                    MaxWS,                                          &
!AQI, 
                    D_TEMP, D_PPFD, SUNT, SHAT,         &
                    SUNF, SUNP, SHAP,ER, NON_DIMGARMA  )
! PURPOSE: Calculate Vegetation emission activity (EA) for each emission
!		class as the product of EA for individual drivers
!		calculated by the following functions
!
! Vegetation Emission Activity (EA) algorithm FUNCTIONS
!
!   GAMTLD: EA Temperature response (light dependent emission) 
!   GAMTLI: EA Temperature response (light independent emission)
!   GAMP: EA Light response
!   GAMTP: combines GAMLD, GAMLI, GAMP to get canopy average
!
!   GAMLA: EA leaf age response 
!   GAMBD: EA bidirectional exchange LAI response
!
!   CDEA: Canopy depth emission response
!
!   GAMHW: EA response to high wind storms
!   GAMHT: EA resposne to high temperature
!   GAMLT: EA response to low temperature
!
!   GAMAQ: EA response to air pollution
!
!   GAMCO2: EA CO2 response (only applied to isoprene)
!   GAMSM: EA response to soil moisture (multiplied with LDF)
!
! INCLUDE FILES
!     'PARMS3.EXT'   ! I/O API parameters
!     'IODECL3.EXT'  ! I/O API function declarations
!     'FDESC3.EXT'   ! I/O API file desc. data structures
!     'MEGVEA.EXT'    ! coefficients
!
!  INPUT Files
!	Single value for each location
!		LDF: Light dependent fraction (for categories other than
!		monoterpene, use constant values from MEGVEA.EXT)
!		AQ:  Air Quality indicator
!	Time series for each location
!		LAI: Leaf Area Index 
!		GAMSM: emission activity response to soil moisture
!		MaxT: Daily maximum temperature (K)
!		MinT: Daily minimum temperature (K)
!		MaxWS: Daily mximum wind speed (m/s)
!               D_TEMP: Daily average temperature (K)
!               D_PPFD: Daily averaged PPFD (umol/m2.s)
!	Hourly time series for each canopy layer in each location
!		sunfrac: fraction of leaves that are sunlit
!		SUNT: leaf Temperature of sunlit leaves (K)
!		SUNP: sun leaf PPFD visible light (micromol/m2/s)
!		SHAT: leaf Temperature of shade leaves (K)
!		SHAP: shade leaf PPFD visible light (micromol/m2/s)
!	
!  OUTPUT Files
!	Hourly time series for each location 
!		Emission activity for each of the 20 emission types
!
!
! HISTORY:
!   Based on code initiated by Alex Guenther in 1990s
!   Coded in FORTRAN as
!   MEGAN: Jack Chen 11/04
!   MEGANv2.04: Tan 11/21/06
!   MEGANv2.1: X. Wang 11/04/2007
!       Modified by Julia Lee-Taylor 03/18/2008
!       Modified by Xuemei Wang 09/30/2008
!       Modified by Tan 07/28/2011
!   MEGAN3.0:
!   Alex Guenther and Ling Huang Feb 2017
!     converted program to subroutine

    IMPLICIT NONE
! INPUT VARIABLES
    INTEGER, INTENT(IN)  ::  LAYERS
    INTEGER, INTENT(IN)  ::  JDATE
    REAL, INTENT(IN)     ::  ZTIME
    REAL, INTENT(IN)     ::  LAIp       ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  LAIc       ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  LDF_in      ( NCOLS, NROWS , 4 ) !only 4 use maps
    REAL, INTENT(IN)     ::  GAMSM_in       ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  MaxT        ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  MinT        ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  MaxWS       ( NCOLS, NROWS )
!    REAL, INTENT(IN)     ::  AQI         ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  D_TEMP      ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  D_PPFD      ( NCOLS, NROWS ) ! comes in as rgrnd
    REAL, INTENT(IN)     ::  SUNT        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SHAT        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SUNF        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SUNP        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SHAP       ( NCOLS, NROWS, LAYERS )

! OUTPUT VARIABLES
    REAL, INTENT(OUT)     :: ER(NCOLS, NROWS )
    REAL, INTENT(OUT)     :: NON_DIMGARMA (NCLASS,NCOLS, NROWS)

 !LOCAL VARIABLES

    LOGICAL, PARAMETER    :: GAMBD_YN  = .false.

    LOGICAL, PARAMETER    :: GAMAQ_YN  = .false.
! For the CMAQ implementation of MEGAN  we refer to soil moisture 
! at layer 2, which is 1 meter for PX and 0.5 m for NOAH.
! Keep this in mind when enabling the GAMSM stress.
    LOGICAL, PARAMETER    :: GAMSM_YN  = .false. 
    LOGICAL, PARAMETER    :: GAMHT_YN  = .false.
    LOGICAL, PARAMETER    :: GAMLT_YN  = .false.
    LOGICAL, PARAMETER    :: GAMHW_YN  = .false.
    LOGICAL, PARAMETER    :: GAMCO2_YN = .false.

    REAL                  :: VPGWT(LAYERS), Ea1L, Ea2L

    REAL  :: CDEA   ( NCOLS, NROWS, LAYERS ) ! Emission response to canopy depth


    REAL  :: YEAR
    REAl  :: DOY
    REAL  :: GAMLA  ( NCOLS, NROWS )     ! EA leaf age response
    REAL  :: GAMAQ  ( NCOLS, NROWS )     ! EA response to air pollution
    REAL  :: GAMBD  ( NCOLS, NROWS )     ! EA bidirectional exchange LAI response
    REAL  :: GAMHT  ( NCOLS, NROWS )     ! EA response to high temperature
    REAL  :: GAMLT  ( NCOLS, NROWS )     ! EA response to low temperature

    REAL  :: GAMHW  ( NCOLS, NROWS )     ! EA response to high wind speed
    REAL  :: GAMSM  ( NCOLS, NROWS )     ! EA response to soil moisture
    REAL  :: GAMCO2 ( NCOLS, NROWS )     ! EA response to CO2
    REAL  :: GAMTP                       ! combines GAMLD, GAMLI, GAMP to get canopy average

    REAL  :: LDFMAP ( NCOLS, NROWS )     ! light depenedent fraction map


    REAL ::  SUM1, SUM2




    ! loop indices
    !INTEGER :: IDATE, ITIME
    INTEGER :: S, T, I, J, K 
!         YEAR   = FLOAT( MOD( JDATE, 1000 ) )
!         DOY    = FLOAT( JDATE / 1000  )       
! these are backwards I think. don't need them anyway

            

    ! EA response to canopy temperature/light

    IF ( Layers .EQ. 5 ) THEN
        VPGWT(1) = 0.1184635
        VPGWT(2) = 0.2393144
        VPGWT(3) = 0.284444444
        VPGWT(4) = 0.2393144
        VPGWT(5) = 0.1184635
    ELSE
        DO K = 1,Layers
            VPGWT(K) = 1.0 / FLOAT( Layers )
        END DO
    ENDIF

! First process Factors independent of species emission classes S :
    
    ! Emission response to canopy depth
    CALL GAMMA_CD( NCOLS, NROWS, Layers, LAIc, CDEA )

    ! EA bidirectional exchange LAI response
    IF ( GAMBD_YN ) THEN
        CALL GAMMA_LAIbidir(NCOLS, NROWS, LAIc, GAMBD)
    ELSE
        GAMBD = 1.0
    ENDIF

    IF ( GAMCO2_YN ) THEN
        CALL GAMMA_CO2(NCOLS, NROWS, GAMCO2)
    ELSE
        GAMCO2 = 1.0
    ENDIF

!  Now process all factors dependent on S:

    DO S = 1,NEMIS  ! Loop over all the emission classes

        IF ( S .EQ. 3 .OR. S .EQ. 4 .OR. S .EQ. 5 .OR. S .EQ. 6 ) THEN
!    otherwise use the input values.           
            LDFMAP = LDF_IN(:,:,S-2) ! only LDF 3, 4, 5, and 6 in file
        ELSE
!  For these species,  Read LDF from previous MEGVEA.EXT 
            LDFMAP = LDF(S)

        ENDIF


        ! leaf age activity factor:  dependent upon S
        CALL GAMMA_A( NCOLS, NROWS, S, LAIp, LAIc, D_TEMP, GAMLA )

        ! emission activity response to air quality

        IF ( GAMAQ_YN ) THEN
!            CALL GAMMA_AQ(NCOLS, NROWS, S, AQI, GAMAQ)
        ELSE
            GAMAQ = 1.0
        ENDIF

        IF ( GAMSM_YN ) THEN
            GAMSM = GAMSM_in
        ELSE
            GAMSM = 1.0
        ENDIF

        ! EA response to high temperature
        IF ( GAMHT_YN ) THEN
            CALL GAMMA_HT(NCOLS, NROWS, S, MaxT, GAMHT)
        ELSE
            GAMHT = 1.0
        ENDIF

        ! EA response to low temperature
        IF ( GAMLT_YN ) THEN
            CALL GAMMA_LT(NCOLS, NROWS, S, MinT, GAMLT)
        ELSE
            GAMLT = 1.0
        ENDIF

        ! EA response to high wind speed
        IF ( GAMHW_YN ) THEN
            CALL GAMMA_HW(NCOLS, NROWS, S, MaxWS, GAMHW)
        ELSE
            GAMHW = 1.0
        ENDIF


!

! THe following code has bee revised for efficiency and clarity.
!          DO I = 1, NCOLS
!            DO J = 1, NROWS
!
!              DO K = 1, Layers
!
!
!                Ea1L(K)  = CDEA(I,J,K) *                              &
!                      GAMTLD(SunT(I,J,K),D_TEMP(I,J),S) *             &
!                      GAMP(SunP(I,J,K),D_PPFD(I,J)) *  SunF(I,J,K) +  &
!                      GAMTLD(ShaT(I,J,K),D_TEMP(I,J),S) *             &
!                      GAMP(ShaP(I,J,K),D_PPFD(I,J))                   &
!                      * (1-SunF(I,J,K))
!
!                Ea2L(K) = GAMTLI(SunT(I,J,K),S)* SunF(I,J,K)+         &
!                   GAMTLI(ShaT(I,J,K),S)*(1-SunF(I,J,K))
!
!              ENDDO   ! ENDDO canopy layers
!              
!            GAMTP(I,J)=SUM((Ea1L(:)*LDFMAP(I,J) +                   &
!                            Ea2L(:)*(1-LDFMAP(I,J)))* VPGWT( : )) )
!            ENDDO   ! NROWS
!          ENDDO ! NCOLS



  
        DO J = 1, NROWS
        DO I = 1, NCOLS! Preserve stride 1 for output arrays

            SUM1 = 0.0
            SUM2 = 0.0

            DO K = 1, Layers
! 2.025 is the conversion to PPFD. 
! SWDNB*.45 = PAR (Wm-2)
! PAR*4.5 = PPFD (umol/m2/s)
                Ea1L = CDEA(I,J,K) *                                  &
                      GAMTLD(SunT(I,J,K),D_TEMP(I,J),S) *             &
              GAMP(SunP(I,J,K),D_PPFD(I,J)*2.025) *  SunF(I,J,K) +    &
                      GAMTLD(ShaT(I,J,K),D_TEMP(I,J),S) *             &
                  GAMP(ShaP(I,J,K),D_PPFD(I,J)*2.025)                 &
                      * (1.0-SunF(I,J,K))
                SUM1 = SUM1 + Ea1L*VPGWT(K)

                Ea2L = GAMTLI(SunT(I,J,K),S)* SunF(I,J,K) +           &
                      GAMTLI(ShaT(I,J,K),S)*(1.0-SunF(I,J,K))
                SUM2 = SUM2 + Ea2L*VPGWT(K)

            END DO   ! END DO canopy layers

            GAMTP = SUM1*LDFMAP(I,J) + SUM2*( 1.0-LDFMAP(I,J) )

 ! ... Calculate emission activity factors
 
            IF ( S .EQ. 1 ) THEN
!    GAMCO2 only applied to isoprene
                ER(:,:) = LAIc(I,J) * GAMTP * GAMCO2(I,J) * GAMLA(I,J) *       &
                          GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) *  &
                          GAMSM(I,J) 

            ELSE IF ( S .EQ. 13 ) THEN
            
 !   GAMBD only applied to ethanol and acetaldehyde
            
                ER(I,J) = LAIc(I,J) * GAMTP * GAMBD(I,J) * GAMLA(I,J) *        &
                      GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) *      &
                      GAMSM(I,J) 

            ELSE
!  Process remaining species            
            
                ER(I,J) = LAIc(I,J) * GAMTP * GAMLA(I,J) *                            &
                      GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) *             &
                       GAMSM(I,J) 

            END IF
            IF ( ER(I,J).GT.0.0 ) THEN
                NON_DIMGARMA (S,I,J) = ER(I,J)
            ELSE
                NON_DIMGARMA (S,I,J) = 0.0
            END IF

        END DO   ! NCOLS
        END DO ! NROWS

    END DO  ! End loop of species (S)
 
    RETURN
    
   END SUBROUTINE MEGVEA


end module megan_fx
