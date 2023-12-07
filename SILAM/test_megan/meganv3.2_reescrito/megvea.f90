module meg_vea
use megan_ini
   implicit none
   !!!Number of emission classes
   !!INTEGER, PARAMETER :: NCLASS = 19
   !!INTEGER, PARAMETER :: NEMIS  = NCLASS
   !!! number of emission classes

   ! CO2 related emission activity factor parameters
   real,parameter :: CO2   = 400.0
   real,parameter :: ISmax =   1.344
   real,parameter :: CO2h  =   1.4614
   real,parameter :: Cstar = 585.0

   !AQI default:
   real,parameter :: AQI_default=10
   ! PSTD
   REAL,PARAMETER :: PSTD = 200

   ! canopy depth emission response
   real,parameter :: CCD1 = -0.2
   real,parameter :: CCD2 =  1.3

   !Light and temperature emission activity response coefficients for each emission class
   real LDF(NCLASS)  !LDF: light dependent fraction
   real CT1(NCLASS)   !CT1: temperature coefficient (emission type 1: light dependent)
   real Cleo(NCLASS)   !Cleo: temperature coefficient (emission type 1: light dependent)
   real beta(NCLASS)   !beta: temperature coefficient (emission type 2: light independent)
      data    beta(1),LDF(1),CT1(1),Cleo(1)        / 0.13,1.0,95,2    /
      data    beta(2),LDF(2),CT1(2),Cleo(2)        / 0.13,1.0,95,2    /
      data    beta(3),LDF(3),CT1(3),Cleo(3)        / 0.10,0.6,80,1.83 /
      data    beta(4),LDF(4),CT1(4),Cleo(4)        / 0.10,0.9,80,1.83 /
      data    beta(5),LDF(5),CT1(5),Cleo(5)        / 0.10,0.2,80,1.83 /
      data    beta(6),LDF(6),CT1(6),Cleo(6)        / 0.10,0.4,80,1.83 /
      data    beta(7),LDF(7),CT1(7),Cleo(7)        / 0.10,0.1,80,1.83 /
      data    beta(8),LDF(8),CT1(8),Cleo(8)        / 0.10,0.0,80,1.83 /
      data    beta(9),LDF(9),CT1(9),Cleo(9)        / 0.17,0.5,130,2.37/
      data    beta(10),LDF(10),CT1(10),Cleo(10)    / 0.17,0.4,130,2.37/
      data    beta(11),LDF(11),CT1(11),Cleo(11)    / 0.08,0.8,60,1.6  /
      data    beta(12),LDF(12),CT1(12),Cleo(12)    / 0.10,0.2,80,1.83 /
      data    beta(13),LDF(13),CT1(13),Cleo(13)    / 0.13,0.8,95,2    /
      data    beta(14),LDF(14),CT1(14),Cleo(14)    / 0.13,0.8,95,2    /
      data    beta(15),LDF(15),CT1(15),Cleo(15)    / 0.10,0.2,80,1.83 /
      data    beta(16),LDF(16),CT1(16),Cleo(16)    / 0.10,0.2,80,1.83 /
      data    beta(17),LDF(17),CT1(17),Cleo(17)    / 0.10,0.8,80,1.83 /
      data    beta(18),LDF(18),CT1(18),Cleo(18)    / 0.10,0.1,80,1.83 /
      data    beta(19),LDF(19),CT1(19),Cleo(19)    / 0.08,1.0,60,1.6  /

   ! Parameters for leaf age algorithm for each emission activity classes
   real Anew(NCLASS)
   real Agro(NCLASS)
   real Amat(NCLASS)
   real Aold(NCLASS)
      data  Anew( 1), Agro( 1), Amat( 1), Aold( 1) / 0.05 ,0.6, 1.0, 0.9  /
      data  Anew( 2), Agro( 2), Amat( 2), Aold( 2) / 0.05 ,0.6, 1.0, 0.9  /
      data  Anew( 3), Agro( 3), Amat( 3), Aold( 3) / 2.0  ,1.8, 1.0, 1.05 /
      data  Anew( 4), Agro( 4), Amat( 4), Aold( 4) / 2.0  ,1.8, 1.0, 1.05 /
      data  Anew( 5), Agro( 5), Amat( 5), Aold( 5) / 2.0  ,1.8, 1.0, 1.05 /
      data  Anew( 6), Agro( 6), Amat( 6), Aold( 6) / 2.0  ,1.8, 1.0, 1.05 /
      data  Anew( 7), Agro( 7), Amat( 7), Aold( 7) / 2.0  ,1.8, 1.0, 1.05 /
      data  Anew( 8), Agro( 8), Amat( 8), Aold( 8) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew( 9), Agro( 9), Amat( 9), Aold( 9) / 0.4  ,0.6, 1.0, 0.95 /
      data  Anew(10), Agro(10), Amat(10), Aold(10) / 0.4  ,0.6, 1.0, 0.95 /
      data  Anew(11), Agro(11), Amat(11), Aold(11) / 3.5  ,3.0, 1.0, 1.2  /
      data  Anew(12), Agro(12), Amat(12), Aold(12) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(13), Agro(13), Amat(13), Aold(13) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(14), Agro(14), Amat(14), Aold(14) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(15), Agro(15), Amat(15), Aold(15) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(16), Agro(16), Amat(16), Aold(16) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(17), Agro(17), Amat(17), Aold(17) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(18), Agro(18), Amat(18), Aold(18) / 1.0  ,1.0, 1.0, 1.0  /
      data  Anew(19), Agro(19), Amat(19), Aold(19) / 1.0  ,1.0, 1.0, 1.0  /

   !stress emission activity response coefficients for each emission class
   REAL    CAQ(NCLASS)   !CAQ: coefficient for poor Air Quality stress
   REAL    CHW(NCLASS)   !CHW: coefficient for high wind speed stress
   REAL    CHT(NCLASS)   !CHT: coefficient for high temperature stress
   REAL    CLT(NCLASS)   !CLT: coefficient for high temperature stress
      DATA  CAQ(1) ,CHW(1) ,CHT(1) ,CLT(1)    / 1,1,1,1 /
      DATA  CAQ(2) ,CHW(2) ,CHT(2) ,CLT(2)    / 1,1,1,1 /
      DATA  CAQ(3) ,CHW(3) ,CHT(3) ,CLT(3)    / 1,5,1,1 /
      DATA  CAQ(4) ,CHW(4) ,CHT(4) ,CLT(4)    / 5,5,5,5 /
      DATA  CAQ(5) ,CHW(5) ,CHT(5) ,CLT(5)    / 1,5,1,1 /
      DATA  CAQ(6) ,CHW(6) ,CHT(6) ,CLT(6)    / 1,5,1,1 /
      DATA  CAQ(7) ,CHW(7) ,CHT(7) ,CLT(7)    / 1,5,1,1 /
      DATA  CAQ(8) ,CHW(8) ,CHT(8) ,CLT(8)    / 1,1,1,1 /
      DATA  CAQ(9) ,CHW(9) ,CHT(9) ,CLT(9)    / 5,5,5,5 /
      DATA  CAQ(10),CHW(10),CHT(10),CLT(10)   / 5,5,5,5 /
      DATA  CAQ(11),CHW(11),CHT(11),CLT(11)   / 1,1,1,1 /
      DATA  CAQ(12),CHW(12),CHT(12),CLT(12)   / 1,1,1,1 /
      DATA  CAQ(13),CHW(13),CHT(13),CLT(13)   / 1,1,1,1 /
      DATA  CAQ(14),CHW(14),CHT(14),CLT(14)   / 1,1,1,1 /
      DATA  CAQ(15),CHW(15),CHT(15),CLT(15)   / 1,1,1,1 /
      DATA  CAQ(16),CHW(16),CHT(16),CLT(16)   / 1,1,1,1 /
      DATA  CAQ(17),CHW(17),CHT(17),CLT(17)   / 5,5,5,5 /
      DATA  CAQ(18),CHW(18),CHT(18),CLT(18)   / 1,1,1,1 /
      DATA  CAQ(19),CHW(19),CHT(19),CLT(19)   / 1,1,1,1 /

   REAL    TAQ(NCLASS)    !TAQ: threshold for poor Air Quality stress (ppm-hours)
   REAL    THW(NCLASS)   !THW: threshold for high wind speed stress (m/s)
   REAL    THT(NCLASS)   !THT: threshold for high temperature stress (Celsius degree)
   REAL    TLT(NCLASS)   !TLT: threshold for high temperature stress (Celsius degree)
      DATA    TAQ(1),THW(1),THT(1),TLT(1)           / 20,12,40,10 /
      DATA    TAQ(2),THW(2),THT(2),TLT(2)           / 20,12,40,10 /
      DATA    TAQ(3),THW(3),THT(3),TLT(3)           / 20,12,40,10 /
      DATA    TAQ(4),THW(4),THT(4),TLT(4)           / 20,12,40,10 /
      DATA    TAQ(5),THW(5),THT(5),TLT(5)           / 20,12,40,10 /
      DATA    TAQ(6),THW(6),THT(6),TLT(6)           / 20,12,40,10 /
      DATA    TAQ(7),THW(7),THT(7),TLT(7)           / 20,12,40,10 /
      DATA    TAQ(8),THW(8),THT(8),TLT(8)           / 20,12,40,10 /
      DATA    TAQ(9),THW(9),THT(9),TLT(9)           / 20,12,40,10 /
      DATA    TAQ(10),THW(10),THT(10),TLT(10)       / 20,12,40,10 /
      DATA    TAQ(11),THW(11),THT(11),TLT(11)       / 20,12,40,10 /
      DATA    TAQ(12),THW(12),THT(12),TLT(12)       / 20,12,40,10 /
      DATA    TAQ(13),THW(13),THT(13),TLT(13)       / 20,12,40,10 /
      DATA    TAQ(14),THW(14),THT(14),TLT(14)       / 20,12,40,10 /
      DATA    TAQ(15),THW(15),THT(15),TLT(15)       / 20,12,40,10 /
      DATA    TAQ(16),THW(16),THT(16),TLT(16)       / 20,12,40,10 /
      DATA    TAQ(17),THW(17),THT(17),TLT(17)       / 20,12,40,10 /
      DATA    TAQ(18),THW(18),THT(18),TLT(18)       / 20,12,40,10 /
      DATA    TAQ(19),THW(19),THT(19),TLT(19)       / 20,12,40,10 /

    !stress emission activity delta thresholds for each emission class
    !DTAQ: delta threshold for poor Air Quality stress (ppm-hours)
    REAL           DTAQ(NCLASS)
    !DTHW: delta threshold for high wind speed stress (m/s)
    REAL           DTHW(NCLASS)
    !DTHT: delta threshold for high temperature stress (Celsius degree)
    REAL           DTHT(NCLASS)
    !DTLT: delta threshold for low temperature stress (Celsius degree)
    REAL           DTLT(NCLASS)
      DATA  DTAQ(1) ,DTHW(1) ,DTHT(1) ,DTLT(1)   / 30,8,8,8 /
      DATA  DTAQ(2) ,DTHW(2) ,DTHT(2) ,DTLT(2)   / 30,8,8,8 /
      DATA  DTAQ(3) ,DTHW(3) ,DTHT(3) ,DTLT(3)   / 30,8,8,8 /
      DATA  DTAQ(4) ,DTHW(4) ,DTHT(4) ,DTLT(4)   / 30,8,8,8 /
      DATA  DTAQ(5) ,DTHW(5) ,DTHT(5) ,DTLT(5)   / 30,8,8,8 /
      DATA  DTAQ(6) ,DTHW(6) ,DTHT(6) ,DTLT(6)   / 30,8,8,8 /
      DATA  DTAQ(7) ,DTHW(7) ,DTHT(7) ,DTLT(7)   / 30,8,8,8 /
      DATA  DTAQ(8) ,DTHW(8) ,DTHT(8) ,DTLT(8)   / 30,8,8,8 /
      DATA  DTAQ(9) ,DTHW(9) ,DTHT(9) ,DTLT(9)   / 30,8,8,8 /
      DATA  DTAQ(10),DTHW(10),DTHT(10),DTLT(10)  / 30,8,8,8 /
      DATA  DTAQ(11),DTHW(11),DTHT(11),DTLT(11)  / 30,8,8,8 /
      DATA  DTAQ(12),DTHW(12),DTHT(12),DTLT(12)  / 30,8,8,8 /
      DATA  DTAQ(13),DTHW(13),DTHT(13),DTLT(13)  / 30,8,8,8 /
      DATA  DTAQ(14),DTHW(14),DTHT(14),DTLT(14)  / 30,8,8,8 /
      DATA  DTAQ(15),DTHW(15),DTHT(15),DTLT(15)  / 30,8,8,8 /
      DATA  DTAQ(16),DTHW(16),DTHT(16),DTLT(16)  / 30,8,8,8 /
      DATA  DTAQ(17),DTHW(17),DTHT(17),DTLT(17)  / 30,8,8,8 /
      DATA  DTAQ(18),DTHW(18),DTHT(18),DTLT(18)  / 30,8,8,8 /
      DATA  DTAQ(19),DTHW(19),DTHT(19),DTLT(19)  / 30,8,8,8 /

contains

subroutine megvea(  ncols,nrows,layers,          &
                    laip, laic,ldf_in,           &
                    GAMSM_in,                    &
                    MaxT, MinT, MaxWS,           &
                    D_TEMP, D_PPFD,              &
                    SUNT, SHAT, SUNF, SUNP, SHAP,&
                    ER, NON_DIMGARMA             )!AQI, 
    implicit none
    ! input variables
    integer, intent(in)  ::  ncols,nrows,layers
    real, intent(in)     ::  laip(ncols,nrows),laic(ncols,nrows)
    real, intent(in)     ::  ldf_in(ncols,nrows,4) !only 4 use maps
    REAL, INTENT(IN)     ::  GAMSM_in(NCOLS,NROWS)
    REAL, INTENT(IN)     ::  MaxT(NCOLS,NROWS),MinT(NCOLS,NROWS),MaxWS(NCOLS,NROWS)
    REAL, INTENT(IN)     ::  D_TEMP(NCOLS,NROWS)
    REAL, INTENT(IN)     ::  D_PPFD(NCOLS,NROWS)      ! comes in as rgrnd
    real, intent(in)     ::  sunt(ncols,nrows,layers)
    real, intent(in)     ::  shat(ncols,nrows,layers)
    real, intent(in)     ::  sunf(ncols,nrows,layers)
    real, intent(in)     ::  sunp(ncols,nrows,layers)
    real, intent(in)     ::  shap(ncols,nrows,layers)
    !REAL, INTENT(IN)     ::  AQI         ( NCOLS, NROWS )

    ! output variables
    real, intent(out)     :: ER(ncols,nrows)                    !emission rate
    real, intent(out)     :: non_dimgarma (ncols,nrows,nclass)  !

    ! local variables
    LOGICAL, PARAMETER    :: GAMBD_YN  = .true. !.false.
    LOGICAL, PARAMETER    :: GAMAQ_YN  = .true. !.false.
    ! For the CMAQ implementation of MEGAN  we refer to soil moisture 
    ! at layer 2, which is 1 meter for PX and 0.5 m for NOAH.
    ! Keep this in mind when enabling the GAMSM stress.
    LOGICAL, PARAMETER    :: GAMSM_YN  = .true. !.false. 
    LOGICAL, PARAMETER    :: GAMHT_YN  = .true. !.false.
    LOGICAL, PARAMETER    :: GAMLT_YN  = .true. !.false.
    LOGICAL, PARAMETER    :: GAMHW_YN  = .true. !.false.
    LOGICAL, PARAMETER    :: GAMCO2_YN = .true. !.false.


    REAL  :: CDEA(LAYERS)  ! Emission response to canopy depth
    REAL  :: GAMLA      ! EA leaf age response
    REAL  :: GAMAQ      ! EA response to air pollution
    REAL  :: GAMBD      ! EA bidirectional exchange LAI response
    REAL  :: GAMHT      ! EA response to high temperature
    REAL  :: GAMLT      ! EA response to low temperature
    REAL  :: GAMHW      ! EA response to high wind speed
    REAL  :: GAMSM      ! EA response to soil moisture
    REAL  :: GAMCO2     ! EA response to CO2
    REAL  :: GAMTP      ! combines GAMLD, GAMLI, GAMP to get canopy average
    REAL  :: LDFMAP     ! light depenedent fraction map

    REAL :: VPGWT(LAYERS)
    REAL ::  SUM1, SUM2,Ea1L,Ea2L
    ! loop indices
    !INTEGER :: IDATE, ITIME
    integer :: s, t, i, j, k 
            
    ! EA response to canopy temperature/light
    IF ( Layers .EQ. 5 ) THEN
        VPGWT(1) = 0.1184635
        VPGWT(2) = 0.2393144
        VPGWT(3) = 0.284444444
        VPGWT(4) = 0.2393144
        VPGWT(5) = 0.1184635
    ELSE
        do k = 1,layers
            VPGWT(K) = 1.0 / FLOAT( Layers )
        end do
    ENDIF


   do j = 1, NROWS
      do i = 1, NCOLS! preserve stride 1 for output arrays
                                                           
       ! First process Factors independent of species emission classes S :
       
       cdea(:)=gamma_cd(layers,laic(i,j))  ! Emission response to canopy depth

       IF ( GAMBD_YN ) THEN
           GAMBD=GAMMA_LAIbidir(LAIc(i,j)) ! EA bidirectional exchange LAI response
       ELSE
           GAMBD = 1.0
       ENDIF

       IF ( GAMCO2_YN ) THEN               ! EA CO2 response
           GAMCO2=GAMMA_CO2(co2)
       ELSE
           GAMCO2 = 1.0
       ENDIF

       !  Now process all factors dependent on S:
       do s = 1,nemis  ! Loop over all the emission classes

           IF ( S .EQ. 3 .OR. S .EQ. 4 .OR. S .EQ. 5 .OR. S .EQ. 6 ) THEN
               LDFMAP = LDF_IN(i,j,S-2) ! only LDF 3, 4, 5, and 6 in file
           ELSE
               LDFMAP = LDF(S) !For these species,  Read LDF from previous MEGVEA.EXT 
           ENDIF

           ! leaf age activity factor:  dependent upon S
           GAMLA=GAMMA_A(S, LAIp(i,j), LAIc(i,j), D_TEMP(i,j))

           ! EA response to air quality
           IF ( GAMAQ_YN ) THEN
              GAMAQ=GAMMA_AQ(S, AQI_default)
           ELSE
               GAMAQ = 1.0
           ENDIF

           IF ( GAMSM_YN ) THEN
               GAMSM = GAMSM_in(i,j)
           ELSE
               GAMSM = 1.0
           ENDIF
           ! EA response to high temperature
           IF ( GAMHT_YN ) THEN
               GAMHT=GAMMA_HT(S, MaxT(i,j))
           ELSE
               GAMHT = 1.0
           ENDIF
           ! EA response to low temperature
           IF ( GAMLT_YN ) THEN
               GAMLT=GAMMA_LT(S, MinT(i,j))
           ELSE
               GAMLT = 1.0
           ENDIF
           ! EA response to high wind speed
           IF ( GAMHW_YN ) THEN
             GAMHW=GAMMA_HW(S, MaxWS(i,j))
           ELSE
               GAMHW = 1.0
           ENDIF
           SUM1 = 0.0
           SUM2 = 0.0
           do k = 1, layers
              ! 2.025 is the conversion to PPFD. 
              ! SWDNB*.45 = PAR (Wm-2)
              ! PAR*4.5 = PPFD (umol/m2/s)
             Ea1L = CDEA(K) *                                              &
                   GAMTLD(SunT(i,j,k),D_TEMP(i,j),S) *  GAMP(SunP(i,j,k),D_PPFD(i,j)*2.025) *        SunF(i,j,k)  + &
                   GAMTLD(ShaT(i,j,k),D_TEMP(i,j),S) *  GAMP(ShaP(i,j,k),D_PPFD(i,j)*2.025) * (1.0 - SunF(i,j,k) )
             SUM1 = SUM1 + Ea1L*VPGWT(K)

             Ea2L = GAMTLI(SunT(i,j,k),S) *      SunF(i,j,k)     + &
                    GAMTLI(ShaT(i,j,k),S) * (1.0-SunF(i,j,k))
             SUM2 = SUM2 + Ea2L*VPGWT(K)
           end do   ! end do canopy layers
           GAMTP = SUM1*LDFMAP + SUM2*( 1.0-LDFMAP )
           ! ... Calculate emission activity factors
           ER(I,J) = LAIc(I,J) * GAMTP * GAMLA * GAMHW * GAMAQ* GAMHT * GAMLT *  GAMSM
           IF ( S .EQ. 1 ) THEN
               ER(i,j) =ER(i,j) * GAMCO2  ! GAMCO2 only applied to isoprene
           ELSE IF ( S .EQ. 13 ) THEN   
               ER(i,j) = ER(i,j) * GAMBD  ! GAMBD only applied to ethanol and acetaldehyde
           END IF

           IF ( ER(I,J) .GT. 0.0 ) THEN
               NON_DIMGARMA (i,j,s) = ER(i,j)
           ELSE                   
               NON_DIMGARMA (i,j,s) = 0.0
           END IF
        end do  ! End loop of species (S)
     end do   ! NCOLS
  end do ! NROWS
 
  RETURN
    
END SUBROUTINE MEGVEA

    function gamma_cd(Layers,LAI)  result(cdea)
      implicit none
      integer, intent(in)     :: layers
      real                    :: lai 
      real, dimension(layers) :: cdea
      integer :: k
      do k=1,layers
             cdea(k)=ccd1*min(lai*(k-0.5)/float(layers),3.0)+ccd2
      enddo
      return
    end function
    !----------------------------------------------------------------
    function gamma_laibidir(lai)  result(gambd)
      real, intent(in) :: lai
      real :: gambd

      IF(LAI < 2) THEN
          GAMBD =  0.5 * LAI
      ELSEIF (LAI .LE. 6 ) THEN
          GAMBD = 1 - 0.0625 * (LAI - 2)
      ELSE
          GAMBD = 0.75
      ENDIF
    end function
    !----------------------------------------------------------------
    function gamma_co2(co2)        result(gamco2)

        IMPLICIT NONE
        REAL,INTENT(IN)     :: CO2
        REAL                :: GAMCO2
        ! local
        REAL    :: Ci, cxxx, cyyy

        Ci      = 0.7 * CO2
        IF ( CO2 .EQ. 400.0 ) THEN
            GAMCO2 = 1.0
        ELSE
            cxxx =  Ci**CO2h
            cyyy =  Cstar**CO2h
            GAMCO2 = ISmax- ((ISmax * cxxx) / (cyyy + cxxx))
        END IF

        RETURN
    end function gamma_co2

    !----------------------------------------------------------------
    !
    !   FUNCTION GAMTLD
    !       EA Temperature response (light dependent emission)
    !----------------------------------------------------------------
    FUNCTION GAMTLD(T1,T24,S)
        IMPLICIT NONE
        REAL,PARAMETER :: Ct2 = 230
        INTEGER        :: S
        REAL           :: T1,T24,T240,Topt, X, Eopt, GAMTLD

        T240 = T24

        IF (T1 < 260.0) THEN
            GAMTLD = 0.0
        ELSE
            ! Temperature at which maximum emission occurs
            Topt = 312.5 + 0.6 * (T240 - 297.0)
            X    = ((1.0 / Topt) - (1.0 / T1)) / 0.00831
            ! Maximum emission (relative to emission at 30 C)
            Eopt = Cleo(S) * EXP(0.05 * (T24 - 297.0)) *          &
                  Exp(0.05*(T240-297.0))

            GAMTLD= Eopt * Ct2 * Exp(Ct1(S) * X) /                &
                  (Ct2 - Ct1(S) * (1.0 - EXP(Ct2 * X)))
        ENDIF
    END FUNCTION GAMTLD
    !----------------------------------------------------------------
    !
    !   FUNCTION GAMTLI
    !       EA Temperature response (light independent emission)
    !----------------------------------------------------------------
    function gamtli(temp,s)
        IMPLICIT NONE

        REAL           :: temp, GAMTLI
        REAL,PARAMETER :: Ts = 303.15
        INTEGER        :: S

        GAMTLI = exp( beta(S)*(temp-Ts) )

    end function gamtli
    !----------------------------------------------------------------
    !
    !   FUNCTION GAMP
    !       EA Light response
    !----------------------------------------------------------------
    function gamp(ppfd1,ppfd24)
        implicit none
        real            :: ppfd1, ppfd24, alpha, c1, gamp
        IF (PPFD24 < 0.01) THEN
            GAMP= 0.0
        ELSE
            Alpha  = 0.004
            !        C1     = 0.0468 * EXP(0.0005 * (PPFD24 - PSTD))
            !     &          * (PPFD24 ** 0.6)
            C1 = 1.03
            !        GAMP= (Alpha * C1 * PPFD1) / ((1 + Alpha**2. * PPFD1**2.)**0.5)
        !   use SQRT her for clarity and efficiency
            GAMP= (Alpha * C1 * PPFD1) / SQRT(1.0 + Alpha**2 * PPFD1**2)
        ENDIF
    end function gamp
    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_HT
    !   EA response to high temperature
    !
    !----------------------------------------------------------------
    function GAMMA_HT(S, MaxT)  result(gamht)
        implicit none
        ! input
        integer,intent(in)     :: s
        real,intent(in)        :: maxt
        ! output
        real                   :: gamht
        ! local
        REAL        :: THTK, t1
            THTK = 273.15 + THT(S)
            t1 = THTK + DTHT(S)
            IF (MaxT <= THTK) THEN
                GAMHT = 1.0
            ELSE IF ( MaxT < t1) THEN
                GAMHT = 1.0 + (CHT(S) - 1.0)* (MaxT - THTK)/DTHT(S)
            ELSE
                GAMHT = CHT(S)
            ENDIF
        RETURN
    end function GAMMA_HT
    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_LT
    !   EA response to low temperature
    !
    !----------------------------------------------------------------
    function GAMMA_LT(S, MinT) result(gamlt)
        implicit none
        ! input
        integer,intent(in):: s
        real,intent(in)   :: mint
        ! output
        real              :: gamlt
        ! local
        real         :: tltk, t1

        TLTK = 273.15 + TLT(S)
        t1 = TLTK - DTLT(S)
        IF (MinT >= TLTK) THEN
            GAMLT = 1.0
        ELSE IF ( MinT > t1) THEN
            GAMLT = 1.0 + (CLT(S) - 1.0)* (TLTK - MinT)/DTLT(S)
        ELSE
            GAMLT = CLT(S)
        ENDIF

        RETURN
    end function GAMMA_LT
    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_HW
    !   EA response to high wind speed
    !
    !----------------------------------------------------------------
    function GAMMA_HW(S, MaxWS) result(gamhw)
        implicit none
        ! input
        integer,intent(in) :: S
        real,intent(in)    :: MaxWS
        ! output
        real               :: GAMHW
        ! local
        real        :: t1

            t1 = THW(S) + DTHW(S)
            IF (MaxWS <= THW(S)) THEN
                GAMHW = 1.0
            ELSE IF ( MaxWS < t1) THEN
                GAMHW = 1.0 + (CHW(S) - 1.0)* (MaxWs - THW(S))/ DTHW(S)
            ELSE
                GAMHW = CHW(S)
            ENDIF
        RETURN
    end function gamma_hw

    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_AQ
    !   EA response to air quality
    !
    !----------------------------------------------------------------
    function gamma_aq(s, aqi) result(gamaq)
        implicit none
        ! input
        integer, intent(in)   :: s
        real,    intent(in)   :: aqi
        ! output
        real                :: gamaq
        ! local
        REAL       :: t1
        t1 = TAQ(S) + DTAQ(S)
        IF (AQI <= TAQ(S)) THEN
            GAMAQ = 1.0
        ELSE IF ( AQI < t1) THEN
            GAMAQ = 1.0 + (CAQ(S) - 1.0)* (AQI - TAQ(S))/DTAQ(S)
        ELSE
            GAMAQ = CAQ(S)
        ENDIF

        RETURN
    end function gamma_aq

    function gamma_a(s, laip, laic, tt) result(gamla)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN):: S
        REAL,INTENT(IN)   :: Tt, LAIp, LAIc
        ! output
        REAL              :: GAMLA

        REAL :: Fnew, Fgro, Fmat, Fold
        REAL :: ti,tm
        REAL       :: TSTLEN  
        !Time step of LAI data
        !if (USE_MEGAN_LAI) THEN
        !  TSTLEN = 8.0 ! 8 daily from MEGAN file
        !else
          TSTLEN = 1.0 ! 1 Daily from soilout/metcro
        !end if

        !---------------------------------------------------
        ! local parameter arrays
        !... Calculate foliage fraction
        IF (LAIp .LT. LAIc) THEN

            !        Calculate ti and tm
            IF (Tt .LE. 303.0) THEN
                ti = 5.0 + 0.7*(300-Tt)
            ELSE
                ti = 2.9
            END IF
            tm = 2.3*ti

            !       Calculate Fnew and Fmat, then Fgro and Fold
            !       Fnew
            IF (ti .GE. TSTLEN) THEN
                Fnew = 1.0 - (LAIp/LAIc)
            ELSE
                Fnew = (ti/TSTLEN) * ( 1-(LAIp/LAIc) )
            END IF

            !       Fmat
            IF (tm .GE. TSTLEN) THEN
                Fmat = LAIp/LAIc
            ELSE
                Fmat = (LAIp/LAIc) + ( (TSTLEN-tm)/TSTLEN ) * ( 1-(LAIp/LAIc) )
            END IF
            Fgro = 1.0 - Fnew - Fmat
            Fold = 0.0

        ELSE IF (LAIp .EQ. LAIc) THEN
            Fnew = 0.0
            Fgro = 0.1
            Fmat = 0.8
            Fold = 0.1

        ELSE IF (LAIp .GT. LAIc) THEN
            Fnew = 0.0
            Fgro = 0.0
            Fold = ( LAIp-LAIc ) / LAIp
            Fmat = 1-Fold
        END IF
        !...  Calculate GAMLA
        GAMLA = Fnew*Anew(S) + Fgro*Agro(S) + Fmat*Amat(S) + Fold*Aold(S)

        RETURN
    end function gamma_a

end module meg_vea
