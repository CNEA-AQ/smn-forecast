module megan_sea

contains
    subroutine megsea (IDATE,ITIME,TSTEP,JYEAR,JDAY,     &
                    L_DESID_DIAG,SLTYP, CTF,LAIc, LAT,   &
                    TEMP, SOILM1,SOILM2, SOILT, PRECADJ, &
                    CFNO, CFNOG, GAMSM, GAMNO, BDSNP_NO  )

!***********************************************************************
!   This subroutine computes soil NO emission activity factor and isoprene
!   soil moisture activity using MCIP output variables.
!
!  DESCRIPTION:
!
!     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
!     to estimate NO emissions
!     Information needed to estimate NO emissions
!     Julian Day          (integer)    JDATE
!     Surface Temperature (MCIP field) TA    (K)
!     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (LSOIL)
!          (ratio of volume of water per volume of soil)
!     Soil Temperature    (MCIP field) SOILT (K)         (LSOIL)
!     Soil Type           (MCIP field) ISLTYP            (LSOIL)
!
!     saturation values for soil types (constants)       (LSOIL)
!     FOR PX Version, the Temperature adjustment factor accounts for wet
!     and dry soils
!                and  the precipitation adjustment factor accounts for
!                saturated soils
!     FOR the non-PX version, the basic algorithm remains with a
!     temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse
!     PULSEDATE julian date for the beginning of an NO pulse
!     PULSETIME        time for the beginning of an NO pulse
!
!     The calculation are based on the following paper by J.J. Yienger
!     and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol
!     100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet
!     and dry soils with the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet
!       and dry adjustment is calculated at each grid cell.  A linear 
!       interpolation between the wet and dry adjustment factor is made 
!       using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture 
!       is determined by taking the MCIP soil moisture field and dividing by the
!       saturation value defined for each soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the
!     following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total.
!       THe types of Pulses as described in YL95 were used to estimate
!       the NO emission rate.
!
!    Also see the following paper for more information:
!    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!    by Tom Pierce and Lucille Bender
!
!    REFERENCES
!
!    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and 
!      Nitric Oxide Emissions from Agricultural Proc Proceedings of the Air and Waste 
!      Management Association/U.S. Environmental Protection Agency EMission Inventory 
!      Conference, Raleigh October 26-28, 1999 Raleigh NC
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start, soil type, 
!     SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!
! HISTORY:
!   07/21/11: Imported from SMOKE-BEIS v3.14 for MEGEAN v2.10 (Tan)
!   03/19/17: Make as an indpendent program (MEGSEA) (Ling Huang)
!   03/31/17: Add calculation for soil moisture activity (Ling Huang)
!   06/10/19: Add an option to use BDSNP model to calculate soil NO
!             emissions (Ling Huang) 
!*********************************************************************

      USE BDSNP_MOD
      USE RUNTIME_VARS, ONLY: BDSNP_MEGAN

      IMPLICIT NONE
 
     ! input variables
     INTEGER, INTENT(IN) :: IDATE, ITIME, TSTEP(3)  
     LOGICAL, INTENT( IN ) :: L_DESID_DIAG
     
     INTEGER, INTENT(IN) :: SLTYP  (NCOLS, NROWS)  ! soil type
     REAL, INTENT(IN)    :: JYEAR, JDAY
     REAL, INTENT(IN)    :: CTF( NrTyp, NCOLS, NROWS ) ! Canopy type factor arra
     REAL, INTENT(IN)    :: LAIc( NCOLS, NROWS )    ! Current time step LAI
     REAL, INTENT(IN)    :: LAT (NCOLS, NROWS )    ! Latitude
     REAL, INTENT(IN)    :: TEMP (NCOLS, NROWS)   ! Temperautre (K)

     REAL, INTENT(IN)    :: SOILM1  (NCOLS, NROWS)  ! soil moisture
     REAL, INTENT(IN)    :: SOILM2  (NCOLS, NROWS)  ! soil moisture
     REAL, INTENT(IN)    :: SOILT  (NCOLS, NROWS)  ! soil temperature
     REAL, INTENT(IN)    :: PRECADJ (NCOLS, NROWS)   

     ! output variable
     REAL, INTENT(OUT)   :: CFNO  (NCOLS, NROWS)       ! Emission activity for crop
     REAL, INTENT(OUT)   :: CFNOG  (NCOLS, NROWS)      ! Emission activity for grass
     REAL, INTENT(OUT)   :: GAMSM  (NCOLS, NROWS)      ! Soil moisture activity for isoprene
     REAL, INTENT(OUT)   :: GAMNO  (NCOLS, NROWS)      ! Final NO emission activity
     REAL, INTENT(OUT)   :: BDSNP_NO (NCOLS, NROWS)    ! BDSNP NO emissions(nmol/s/m2)



! Local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name


      INTEGER :: GDAY, GLEN
      INTEGER :: MXLAI,MXCT
      REAL :: t1,wilt,TMO1,TMO2

      LOGICAL :: LSOIL = .TRUE.

      INTEGER :: T,I,J,MM,DD,I_CT
        
      CFNO  = 0.0 ! INITIALIZE
      CFNOG = 0.0 ! INITIALIZE

        IF (BDSNP_MEGAN) THEN

           CALL GET_DATE(JYEAR, JDAY, MM, DD)
           CALL HRNOBDSNP( IDATE,ITIME,TSTEP,MM,            &
                    L_DESID_DIAG,SOILM1,SOILT,SLTYP,LAIc,   &
                                            BDSNP_NO)
        ELSE
           CALL SOILNOX(IDATE,ITIME,NCOLS,NROWS,           &
                     TEMP,LSOIL,SLTYP, SOILM1, SOILT,      &
                     LAIc, LAT, PRECADJ,                   &
                     CFNO, CFNOG )

           DO I = 1,NCOLS
             DO J = 1,NROWS
               CALL GROWSEASON(IDATE,LAT(I,J),GDAY,GLEN)
               IF (GDAY .EQ. 0) THEN
                ! non growing season
                ! CFNOG for everywhere
                  GAMNO(I,J) = CFNOG(I,J)

                ELSE IF (GDAY .GT. 0 .AND. GDAY .LE. 366) THEN
                ! growing season
                ! CFNOG for everywhere except crops
                TMO1 = 0.
                TMO2 = 0.
                DO I_CT = 1,5
                  TMO1 = TMO1 + CTF(I_CT,I,J)
                  TMO2 = TMO2 + CTF(I_CT,I,J) * CFNOG(I,J)
                ENDDO
                ! CFNO for crops
                TMO1 = TMO1 + CTF(6,I,J)
                TMO2 = TMO2 + CTF(6,I,J) * CFNO(I,J)
                IF (TMO1 .EQ. 0.0) THEN
                   GAMNO(I,J) = 0.0
                ELSE
                   GAMNO(I,J) = TMO2 / TMO1
                ENDIF
                ENDIF
 
              ENDDO  !NCOLS
           ENDDO  !NROWS

        END IF ! YL or BDSNP


        !SOIL MOISTURE FOR ISOPRENE:
        DO I = 1, NCOLS
          DO J = 1, NROWS

            !wilt = WWLT(SLTYP(I,J))
            wilt = Grid_Data%WWLT(I,J)
            t1 = wilt + d1
            IF ( SOILM2(I,J) < wilt ) THEN
                GAMSM(I,J) = 0
            ELSE IF ( SOILM2(I,J) >= wilt .AND. SOILM2(I,J) < t1 ) THEN
                GAMSM(I,J) = (SOILM2(I,J) - wilt)/d1
            ELSE
                GAMSM(I,J) = 1
            END IF
          END DO ! NCOLS
        END DO ! NROWS
         
         
  END SUBROUTINE MEGSEA





!=======================================================================
!=======================================================================
      REAL FUNCTION FERTLZ_ADJ( DATE, LAT )

!***********************************************************************
!  DESCRIPTION:
!    This internal function computes a fertilizer adjustment factor
!    for the given date in yyyyddd format. If it is not growing 
!    season, the adjustment factor is 0; otherwise, it ranges from
!    0.0 to 1.0.
!
!  CALL:
!    GROWSEASON
!
!  HISTORY:
!    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
!***********************************************************************

      IMPLICIT NONE
            
!.... Function arguments
      INTEGER, INTENT(IN) :: DATE
      REAL,    INTENT(IN) :: LAT

!.... Local variables
      INTEGER  GDAY, GLEN

      CHARACTER(256)  MESG         ! message buffer
!-----------------------------------------------------------------------------

      CALL GROWSEASON( DATE, LAT, GDAY, GLEN )
          FERTLZ_ADJ = 0. !INITIALIZE
      IF( GDAY == 0 ) THEN
          FERTLZ_ADJ = 0.
      ELSE IF( GDAY >= 1 .AND. GDAY < 30 ) THEN
          ! first month of growing season
          FERTLZ_ADJ = 1.
      ELSE IF( GDAY >= 30 .AND. GDAY <= 366) THEN
          ! later month of growing season
          FERTLZ_ADJ = 1. + 30. / FLOAT(GLEN) - FLOAT(GDAY) / FLOAT(GLEN)
      ENDIF

!******************  FORMAT  STATEMENTS   ******************************
94010 FORMAT( A, F10.2, 1X, A, I3, ',', I3 )


      RETURN

      END FUNCTION FERTLZ_ADJ
!=======================================================================
!=======================================================================


!=======================================================================
!=======================================================================
      REAL FUNCTION VEG_ADJ( LAI )

!***********************************************************************
!  DESCRIPTION
!    This internal function computes a vegetation adjustment factor
!    based on LAIv.  See Yienger and Levy 1995
!    VEG_ADJ = (EXP(-0.24*LAIv)+EXP(-0.0525*LAIv))*0.5 
!
!  CALL
!    NONE
!
!  HISTORY:
!***********************************************************************

      IMPLICIT NONE
      
!...  Function arguments
      REAL,    INTENT(IN) :: LAI

!-----------------------------------------------------------------------------
      VEG_ADJ = 0.0
      VEG_ADJ = (EXP(-0.24*LAI)+EXP(-0.0525*LAI))*0.5 

!******************  FORMAT  STATEMENTS   ******************************

      RETURN
      END FUNCTION VEG_ADJ
!=======================================================================
!=======================================================================
            


!=======================================================================
!=======================================================================

!=======================================================================
!=======================================================================
      SUBROUTINE GROWSEASON ( DATE, LAT, GDAY, GLEN )

!***********************************************************************
!  DESCRIPTION
!    This internal function computes the day of the growing season
!    corresponding to the given date in yyyyddd format.
!
!  CALL
!    G2J
!
!  HISTORY:
!    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
!               Variation of growing season depends on latitude
!               (Guenther)
!    04/22/2019 Converted to 90 format and redid error repotring
!               modified G2j to be completely internal. 
!               DR. FRANCIS S. BINKOWSKI, IE, UNC-CHAPEL HILL
!***********************************************************************

      IMPLICIT NONE

!.......  Function arguments
      INTEGER, INTENT(IN)  :: DATE
      REAL,    INTENT(IN)  :: LAT
      INTEGER, INTENT(OUT) :: GDAY, GLEN 
!.......  External functions

!.......  Local parameters
      INTEGER            :: GSEASON_START
      INTEGER            :: GSEASON_END

!.......  Local variables
      INTEGER  YEAR, MONTH, DAY
      INTEGER  JDAY
      INTEGER  GSJULIAN_START
      INTEGER  GSJULIAN_END
      OPEN(UNIT = 10, FILE = 'growseason_error_outoputs.txt')
!-----------------------------------------------------------------------------
!   NOTE: The use of "julian Day" to describe the day of tHE year is
!     technically incorrect. 

 ! The Julian Day Number (JDN) is the integer assigned to a whole solar 
 ! day in the Julian day count starting from noon Universal time, with 
 ! Julian day number 0 assigned to the day starting at noon on Monday, 
 ! January 1, 4713 BCE, proleptic Julian calendar (November 24, 4714 BCE, 
 ! in the proleptic Gregorian calendar), a date at which three 
 ! multi-year cycles started (which are: Indiction, Solar, and Lunar cycles)
 !  and which preceded any dates in recorded history. 
 !
 !  For example for January 1st, 2000 CE  at 00:00:00.0 UT1 the Julian Day 
 !  is 2451544.500000 according to  the U.S. Naval Observatory.



      YEAR = INT( FLOAT( DATE ) / 1000. )
      JDAY = DATE - YEAR * 1000

      IF ( LAT .LE. 23.0 .AND. LAT .GE. -23.0 ) THEN
      ! tropical regions, year round
         GSEASON_START = 0101
         GSEASON_END   = 1231

         GSJULIAN_START = G2J(YEAR, GSEASON_START)
         GSJULIAN_END   = G2J(YEAR, GSEASON_END)
         GDAY = JDAY - GSJULIAN_START + 1
         GLEN = GSJULIAN_END - GSJULIAN_START + 1
      ELSE IF ( LAT .LT. -23.0 ) THEN
      ! southern hemisphere
         IF ( LAT .LT. -60.0 ) THEN
         ! antarctic start = 0 end = 0, no growing
            GDAY = 0
            GLEN = 0
         ELSE
         ! southern hemisphere temperate, NOV, DEC, JAN-MAY
            IF (JDAY .GE. 1101 .AND. JDAY .LE. 1231 ) THEN
              GSEASON_START = 1101
              GSEASON_END   = 1231

              GSJULIAN_START = G2J(YEAR, GSEASON_START)
              GSJULIAN_END   = G2J(YEAR, GSEASON_END)
              GDAY = JDAY - GSJULIAN_START + 1
            ELSE IF (JDAY .GE. 0101 .AND. JDAY .LE. 0531) THEN
              GSEASON_START = 0101
              GSEASON_END   = 0531

              GSJULIAN_START = G2J(YEAR, GSEASON_START)
              GSJULIAN_END   = G2J(YEAR, GSEASON_END)
              GDAY = JDAY - GSJULIAN_START + 1 + 61
            ELSE
              GDAY = 0
            ENDIF
            GLEN = 30 + 31 + G2J(YEAR,0531) - G2J(YEAR,0101) + 1

         ENDIF
      ELSE IF ( LAT .GT. 23.0 ) THEN
      ! northern hemisphere
         IF ( LAT .GT. 65.0 ) THEN
         ! arctic start = 0 end = 0, no growing season
            GDAY = 0
            GLEN = 0
         ELSE
         ! northern hemisphere temperate
         ! start= (lat-23)*4.5            189
         ! end = 365 -((lat-23)*3.3)      226
            GSEASON_START = 0
            GSEASON_END   = 1231

            GSJULIAN_START = 0
            GSJULIAN_END   = G2J(YEAR, GSEASON_END)

            GSJULIAN_START = INT( (LAT-23.0) * 4.5 )
            GSJULIAN_END   = GSJULIAN_END - INT( (LAT-23.0) * 3.3 )
            
            ! UNC added to avoid GDAY excede 366
            IF ( JDAY == 366 .AND. GSJULIAN_START==0 ) GSJULIAN_END = GSJULIAN_END - 1
                        
            IF (JDAY .GE. GSJULIAN_START .AND. JDAY .LE. GSJULIAN_END) THEN
               GDAY = JDAY - GSJULIAN_START + 1
            ELSE
               GDAY = 0
            ENDIF
            GLEN = GSJULIAN_END - GSJULIAN_START + 1
         ENDIF
      ENDIF
  
     RETURN

      END SUBROUTINE GROWSEASON

!=======================================================================
!=======================================================================
      
!   This is a modified version of function G2J.

      
      INTEGER FUNCTION G2J( YYYY, MMDD )
      IMPLICIT NONE

!.......  Function arguments
      INTEGER, INTENT(IN) :: YYYY
      INTEGER, INTENT(IN) :: MMDD


!.......  Local parameters
      INTEGER :: MM
      INTEGER :: DD
      INTEGER :: K
      INTEGER :: DOY
      LOGICAL :: LEAP


      MM = INT( FLOAT( MMDD ) / 100.0 )
      DD = MMDD - MM * 100

!   use internal code to get G2j

 !     G2J = JULIAN( YYYY, MM , DD )

!   The following code is taken from NASA subroutine get_DOY 
! Original  Programmer:   David G. Simpson
!            NASA Goddard Space Flight Center
!            Greenbelt, Maryland  2077  Date:         November 20, 2001
!  Modified April 13, 2019 by Dr Francis S. Binkowski to do only Gregorian years


          LEAP = .FALSE.
      
!   TEST FOR LEAP YEARS
      
      IF ( MOD(YYYY,4)   .EQ. 0) LEAP = .TRUE.
      IF ( MOD(YYYY,100) .EQ. 0) LEAP = .FALSE.
      IF  (MOD(YYYY,400) .EQ. 0) LEAP = .TRUE.

      IF (LEAP) THEN
         K = 1
      ELSE
         K = 2
      END IF

!   CALCULATE DAY OF THE YEAR

      DOY = ( ( 275 * MM) / 9 ) - K * ( ( MM + 9) / 12 ) + DD - 30
      
      G2J = DOY  

      END FUNCTION G2J

!=======================================================================



      SUBROUTINE SOILNOX( JDATE, JTIME, NX, NY,            &
                      TA, LSOIL, ISLTYP, SOILM, SOILT,     &
                      LAIc, LAT,                           &
                      PRECADJ,                             &
                      CFNO, CFNOG )

!***********************************************************************
!  DESCRIPTION:
!  
!     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
!     to estimate NO emissions 
!     Information needed to estimate NO emissions
!     Julian Day          (integer)    JDATE
!     Surface Temperature (MCIP field) TA    (K)
!     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (LSOIL)
!          (ratio of volume of water per volume of soil)
!     Soil Temperature    (MCIP field) SOILT (K)         (LSOIL)
!     Soil Type           (MCIP field) ISLTYP            (LSOIL)
!
!     saturation values for soil types (constants)       (LSOIL)
!     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
!                and  the precipitation adjustment factor accounts for saturated soils
!     FOR the non-PX version, the basic algorithm remains with a temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse 
!     PULSEDATE julian date for the beginning of an NO pulse 
!     PULSETIME        time for the beginning of an NO pulse
!  
!     The calculation are based on the following paper by J.J. Yienger and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet and dry soils with
!       the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet and dry adjustment is 
!       calculated at each grid cell.  A linear interpolation between the wet and dry adjustment
!       factor is made using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture is determined by
!       taking the MCIP soil moisture field and dividing by the saturation value defined for each
!       soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total. 
!       THe types of Pulses as described in YL95 were used to estimate the NO emission
!       rate.  
!
!    Also see the following paper for more information:
!    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!    by Tom Pierce and Lucille Bender       
!
!    REFERENCES
!
!    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and Nitric Oxide Emissions from Agricultural Processes
!       Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!        Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start,
!     soil type, SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!     
!  REVISION  HISTORY:
!    10/01 : Prototype by GAP
!    10/03 : modified transition to non growing season for jul-oct of the year
!    08/04 : Converted to SMOKE code style by C. Seppanen
!    07/21/11 : Imported form SMOKE-BEIS v3.14 for MEGAN v2.10
!    MAY 13, 2019 made inot f90 format and  improved efficiency - 
! 
!***********************************************************************

!        USE SOILNOX_FX

        IMPLICIT NONE
        

!.........  ARGUMENTS and their descriptions
        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      !  no. columns
        INTEGER, INTENT (IN)  :: NY      !  no. rows

        REAL, INTENT (IN)  ::  TA      ( NX, NY )    !  air temperature (K)
        REAL, INTENT (IN)  ::  SOILM   ( NX, NY )    !  soil moisture (m3/m3)
        REAL, INTENT (IN)  ::  SOILT   ( NX, NY )    !  soil temperature (K)
        REAL, INTENT (IN)  ::  PRECADJ ( NX, NY )    !  precip adjustment
        REAL, INTENT (IN)  ::  LAIc    ( NX, NY )    !  soil temperature (K)
        REAL, INTENT (IN)  ::  LAT     ( NX, NY )    !  Latitude
        REAL, INTENT (IN OUT)  ::  CFNO    ( NX, NY )    !  NO correction factor
        REAL, INTENT (IN OUT)  ::  CFNOG   ( NX, NY )    !  NO correction factor for grass
        
        INTEGER, INTENT (IN)  ::  ISLTYP  ( NX, NY )    !  soil type

        LOGICAL, INTENT (IN) :: LSOIL              ! true: using PX version of MCIP
        
!.........  Local ARRAYS
! Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
!       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
! See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
        INTEGER, PARAMETER :: MAXSTYPES = 16
!        REAL, PARAMETER    :: SATURATION( MAXSTYPES )     =  (/   &       
!                              0.395, 0.410, 0.435, 0.485,         &
!                              0.451, 0.420, 0.477, 0.476,         &
!                              0.426, 0.482, 0.482            /)       

!.........  SCRATCH LOCAL VARIABLES and their descriptions:
        INTEGER       ::   R, C, L      ! counters
        INTEGER       ::   SOILCAT      ! soil category
        
        REAL          ::   CF           ! NO correction factor
        REAL          ::   CFG          ! NO correction factor for grasslands
        REAL          ::  TAIR         ! surface temperature
        REAL          ::   TSOI         ! soil temperature
        REAL          ::   CFNOWET, CFNODRY, RATIO, FAC1, FAC2 
        REAL, PARAMETER ::  const1 = (1. / 3.0 )  * (1.0 / 30.0)
        REAL, PARAMETER ::  const2 =EXP(-0.103 * 30.0)
        CHARACTER(256)  MESG         ! message buffer
        
        CHARACTER(16) :: PROGNAME = 'SOILNOX'   !  program name

!***********************************************************************

 
!.....  Loop through cells
        DO R = 1, NY
        DO C = 1, NX

          TAIR = TA( C, R )         ! unit in degree K

!.......  Check max bounds for temperature

          IF (TAIR > 315.0 ) THEN
              TAIR = 315.0
          END IF

!.......  CFNOG
          IF( TAIR > 303.00 ) TAIR = 303.00

          IF ( TAIR > 268.8690 ) THEN  
              CFG = EXP( 0.04686 * TAIR - 14.30579 ) ! grass (from BEIS2)
          ELSE
              CFG = 0.0
          END IF

          CFNOG(C,R) = CFG
          
!   pre calculate common factors


              FAC2 = const2

!.......  CFNO
          IF( .NOT. LSOIL ) THEN
          ! no soil

             TSOI = 0.72 * TAIR + 82.28
             IF (TSOI <= 273.16) TSOI = 273.16
             IF (TSOI >= 303.16) TSOI = 303.16
             
              FAC1 = (TSOI- 273.16)

             
              
!             CFNODRY = (1./3.) * (1./30.) * (TSOI-273.16)  ! see YL 1995 Equa 9a p. 11452             
              CFNODRY = const1 * FAC1  ! see YL 1995 Equa 9a p. 11452
            
             IF (TSOI <= 283.16) THEN         ! linear cold case
!                 CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
                 CFNOWET =  FAC1 * FAC2 * 0.28 ! see YL 1995 Equ 7b
                 
             ELSE                             ! exponential case
!                 CFNOWET = EXP(0.103 * (TSOI-273.16)) *  EXP(-0.103 * 30.0)
                 CFNOWET = EXP(0.103 * FAC1) *  FAC2
             END IF
             CF = 0.5 * CFNOWET + 0.5 * CFNODRY

          ELSE
          ! soil

             TSOI = SOILT( C,R )
             IF (TSOI <= 273.16) TSOI = 273.16
             IF (TSOI >= 303.16) TSOI = 303.16

              FAC1 = (TSOI- 273.16)

!             CFNODRY = (1./3.)*(1./30.)*(TSOI-273.16)  ! see YL 1995 Equa 9a p. 11452
             CFNODRY = const1 * FAC1  ! see YL 1995 Equa 9a p. 11452
                          
             IF (TSOI <= 283.16) THEN         ! linear cold case
!                CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
                CFNOWET = FAC1 * FAC2 * 0.28 ! see YL 1995 Equ 7b
                
             ELSE                             ! exponential case
!                CFNOWET = EXP(0.103 * (TSOI-273.16)) * EXP(-0.103 * 30.0)
                 CFNOWET = EXP(0.103 * FAC1 ) * FAC2
               
             END IF

             SOILCAT = INT( ISLTYP( C,R ) )
             IF( SOILCAT > 0 .AND. SOILCAT <= MAXSTYPES ) THEN
                 IF(Grid_Data%WSAT(C,R) .eq. 0) THEN
                  ! first ldesid diag call. Do nothing.
                  CF = 0. 
                 ELSE 
                  RATIO = SOILM( C,R ) / Grid_Data%WSAT( C,R )
                  CF = RATIO*CFNOWET + (1.0 - RATIO ) * CFNODRY
                 END IF
             ELSE
             
                 CF = 0.0
                 
             END IF

          END IF  ! Endif LSOIL


!          CFNO(C,R) = CF *                                      &
!                     FERTLZ_ADJ( JDATE, LAT(C,R) ) *           &
!                     VEG_ADJ( LAIc(C,R) ) * PRECADJ(C,R)

          CFNO(C,R) = CF *                                     &
                     FERTLZ_ADJ( JDATE, LAT(C,R) ) *           &
                     VEG_ADJ( LAIc(C,R) ) * PRECADJ(C,R)

          if(cfno(c,r) .lt. 0) then
             cfno(c,r) = 0
          end if

        END DO  ! loop over columns
        END DO  ! loop over rows

        RETURN

        END SUBROUTINE SOILNOX



end module
