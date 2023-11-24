module megcanopy
   implicit none
     real , parameter :: convertwm2toumolm2s = 4.5, solarconstant = 1367., waterairratio = 18.016/28.97
     real,parameter :: Sb = 0.0000000567 ! Stefan-boltzman constant  W m-2 K-4
     ! Canopy characteristics for MEGCAN canopy types
     INTEGER , PARAMETER :: NRTYP = 6          ! Number of canopy types
     INTEGER , PARAMETER :: NRCHA = 17         ! Number of canopy characteristics
     ! 16 variables are assigned for each canopy type
     ! 1  = canopy depth
     ! 2  = leaf width
     ! 3  = leaf length
     ! 4  = canopy height
     ! 5  = scattering coefficient for PPFD
     ! 6  = scattering coefficient for near IR
     ! 7  = reflection coefficient for diffuse PPFD
     ! 8  = reflection coefficient for diffuse near IR
     ! 9  = clustering coefficient (accounts for leaf clumping influence on mean projected leaf area in the direction of the suns beam)
     ! 10 = leaf IR emissivity
     ! 11 = leaf stomata and cuticle factor: 1=hypostomatous, 2=amphistomatous, 1.25=hypostomatous but with some transpiration through cuticle
     ! 12 = daytime temperature lapse rate (K m-1)
     ! 13 = nighttime temperature lapse rate (K m-1)
     ! 14 = warm (>283K) canopy total humidity change (Pa)
     ! 15 = cool (>= 283K) canopy total humidity change (Pa)
     ! 16 = normalized canopy depth where wind is negligible
     ! 17 = canopy transparency
     REAL , DIMENSION(NRCHA,NRTYP) :: Canopychar = reshape( &                                                ! canopy types 
     [16.0, 0.005, 0.10, 24.0, 0.2, 0.8, 0.057, 0.389, 0.85, 0.95, 1.25, 0.06, -0.06, 700. , 150., 0.7, 0.2,&! 1 = Needleleaf trees                      
      16.0, 0.050, 0.10, 24.0, 0.2, 0.8, 0.057, 0.389, 1.10, 0.95, 1.25, 0.06, -0.06, 700. , 150., 0.7, 0.2,&! 2 = Tropical forest trees,
      16.0, 0.050, 0.10, 24.0, 0.2, 0.8, 0.057, 0.389, 0.90, 0.95, 1.25, 0.06, -0.06, 700. , 150., 0.7, 0.2,&! 3 = Temperate broadleaf trees
       1.0, 0.015, 0.10,  2.0, 0.2, 0.8, 0.057, 0.389, 0.85, 0.95, 1.00, 0.06, -0.06, 700. , 150., 0.7, 0.2,&! 4 = shrubs
       0.5, 0.010, 0.15,  0.5, 0.2, 0.8, 0.057, 0.389, 0.70, 0.95, 1.25, 0.06, -0.06, 700. , 150., 0.7, 0.2,&! 5 = herbaceous
       1.0, 0.020, 0.15,  1.0, 0.2, 0.8, 0.057, 0.389, 0.65, 0.95, 1.25, 0.06, -0.06, 700. , 150., 0.7, 0.2]&! 6 = crops
      ,shape=[NRCHA,NRTYP])

contains

SUBROUTINE MEGCAN(YEAR,NCOLS,NROWS,LAYERS,DOY,ZTIME,                &
                 LAT, LONG, LAIc, TEMP, PAR, WIND, PRES, QV, CTF,   &
                 ShadeleafTK, SunleafTK, SunFrac, SunPPFD, ShadePPFD)
!*****************************************************************
! INPUTs
!   Day                  Julian day
!   Lat                  Latitude
!   Long                 Longitude
!   Hour                 Hour of the day
!   TEMP                 Temperature [K]
!   PPFD                 Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind                 Wind speed [m s-1]
!   Humidity             Relative humidity [%]
!   Cantype              Defines set of canopy characteristics
!   LAI                  Leaf area index [m2 per m2 ground area]
!   Pres                 Pressure [Pa]
!*****************************************************************
! OUTPUTs
! For each time step and location
! Each variable is an array with a value for each canopy layer
!                              (vertical profile)
! i = 1 is the top canopy layer, 2 is the next layer, etc.
!   ShadeleafTK          leaf temperature for shade leaves [K] (weighted by canopy type)
!   SunleafTK            leaf temperature for sun leaves [K] (weighted by canopy type)
!   SunFrac              fraction of sun leaves (weighted by canopy type)
!   SunPPFD              PPFD on a sun leaf [umol/m2/s] (weighted by canopy type)
!   ShadePPFD            PPFD on a shade leaf [umol/m2/s] (weighted by canopy type)
     
      IMPLICIT NONE
!   INPUT VARIABLES
       INTEGER, INTENT(IN):: LAYERS,NCOLS,NROWS
       REAL, INTENT(IN)   :: YEAR
       REAl, INTENT(IN)   :: DOY
       REAL, INTENT(IN)   :: ZTIME
       REAL, INTENT(IN)   :: LAT ( NCOLS, NROWS )
       REAL, INTENT(IN)   :: LONG( NCOLS, NROWS )
       REAL, INTENT(IN)   :: LAIc( NCOLS, NROWS )
       REAL, INTENT(IN)   :: TEMP( NCOLS, NROWS )
       REAL, INTENT(IN)   :: PAR ( NCOLS, NROWS ) ! comes in as rgrnd from ASX
       REAL, INTENT(IN)   :: WIND( NCOLS, NROWS )
       REAL, INTENT(IN)   :: PRES( NCOLS, NROWS )
       REAL, INTENT(IN)   :: QV  ( NCOLS, NROWS )
       REAL, INTENT(IN)   :: CTF(NRTYP, NCOLS, NROWS ) ! Canopy type factor array
!   OUTPUT VARIABLES 
       REAL, INTENT(OUT)  ::  ShadeleafTK( NCOLS, NROWS, LAYERS )
       REAL, INTENT(OUT)  ::  SunleafTK  ( NCOLS, NROWS, LAYERS )
       REAL, INTENT(OUT)  ::  SunFrac    ( NCOLS, NROWS, LAYERS )
       REAL, INTENT(OUT)  ::  SunPPFD    ( NCOLS, NROWS, LAYERS )
       REAL, INTENT(OUT)  ::  ShadePPFD  ( NCOLS, NROWS, LAYERS )
!   LOCAL VARIABLES
      INTEGER             :: I, I_CT, J, MM, DD
      INTEGER             :: IDAY         ! For using original solar method

      REAL                :: PPFD!!#(NCOLS, NROWS)
      REAL                :: TotalCT
      REAL                :: month, Date, JDAY
      REAL                :: Sinbeta, Beta, HOUR, DAY
      REAL,DIMENSION(LAYERS) ::  VPgausWt, VPgausDis2,             &
       VPgausDis, VPslwWT,  QdAbsV, QsAbsV, QdAbsn,                &
       QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn,                     &
       TairK, HumidairPa, Ws, SunleafSH, sun_ppfd,shade_ppfd,      &
       SunleafLH,SunleafIR, ShadeleafSH, sun_tk,shade_tk,sun_frac, &
       ShadeleafLH,ShadeleafIR, sun_ppfd_total, shade_ppfd_total,  &
       sun_tk_total, shade_tk_total, sun_frac_total
      REAL :: Solar, Maxsolar, Eccentricity,                       &
              Difffrac, PPFDfrac, QbAbsn,                          &
               Trate, Qbeamv,Qdiffv, Qbeamn, Qdiffn,               &
               QbAbsV,Ea1tCanopy, Ea1pCanopy,                      &
               TairK0, HumidairPa0, Ws0, SH
!      REAL ::  CalcEccentricity,WaterVapPres,                      &      
!               Stability, Calcbeta
   
   
   !loop on each grid-cell:
   do j=1,NCOLS
      do i=1, NROWS
   
      !initialize all vars:
      
      !default values for output
      SunleafTK(I,J,:)   = TEMP(I,J)
      ShadeleafTK(I,J,:) = TEMP(I,J)
      SunPPFD(I,J,:)     = PPFD(I,J)
      ShadePPFD(I,J,:)   = PPFD(I,J)
      SunFrac(I,J,:)     = 1.0
      TotalCT           = 0.0


      TotalCT=sum(CTF(i,j,:))*0.01
      !DO I_CT = 1,NRTYP   !canopy types
      !  TotalCT = TotalCT + CTF(I_CT,I,J) * 0.01
      !ENDDO   ! ENDDO I_CT
   
      if (totalCT .gt. 0 .AND. LAIc(i,j) .gt. 0) then
   
             !(1) calc solar angle
              TairK0   = TEMP(I,J)       !temp (from meteo)   
              Ws0      = WIND(I,J)       !wind (from meteo)   
              Solar    = PPFD(I,J)/2.25  !rad. (from meteo) 

              !	(could it be obtained from COSZEN variable on wrfout)
              Beta     = Calcbeta(Day , Lat(I,J) , Hour )
              Sinbeta  = SIN(Beta  / 57.29578)
              Maxsolar = Sinbeta  * SolarConstant * CalcEccentricity(Day)
   
             !(2) gaussian dist.
             !
             Call GaussianDist(VPgausDis, Layers)
   
             !(3) determine fraction of diffuse PPFD, direct PPFD, diffuse near IR, direct near IR
             !
             Call SolarFractions(Solar,Maxsolar, Qdiffv,Qbeamv,Qdiffn,Qbeamn)
   
             !(4) canopy radiation dist (ppdf)
             !
             Call CanopyRad(VPgausDis, Layers, LAIc(I,J), Sinbeta, Qbeamv, &
                   Qdiffv, Qbeamn, Qdiffn,I_CT ,Canopychar, sun_frac,      &
                   QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv,  &
                   ShadeQv, SunQn, ShadeQn, sun_ppfd, shade_ppfd,          &
                   NrCha,NrTyp)
   
             HumidairPa0  =  WaterVapPres(QV(I,J), PRES(I,J), WaterAirRatio)
             Trate        =  Stability(Canopychar, I_CT, Solar , NrCha, NrTyp)
             !(5) canopy energy balance (temp)
             !
             Call CanopyEB(Trate, Layers, VPgausDis, Canopychar, I_CT,    &
                   TairK, HumidairPa, Ws, sun_ppfd,                       &
                   shade_ppfd, SunQv, ShadeQv, SunQn, ShadeQn,            &
                   sun_tk, SunleafSH, SunleafLH, SunleafIR,               &
                   shade_tk,ShadeleafSH,ShadeleafLH,ShadeleafIR,          &
                   NrCha, NrTyp, Ws0, TairK0, HumidairPa0)
             !(6) compute output variables
             !
             sun_ppfd_total(:)   = sun_ppfd_total(:)   + 0.01*CTF(I_CT,I,J)*sun_ppfd(:)
             shade_ppfd_total(:) = shade_ppfd_total(:) + 0.01*CTF(I_CT,I,J)*shade_ppfd(:)
             sun_tk_total(:)     = sun_tk_total(:)     + 0.01*CTF(I_CT,I,J)*sun_tk(:)
             shade_tk_total(:)   = shade_tk_total(:)   + 0.01*CTF(I_CT,I,J)*shade_tk(:)
             sun_frac_total(:)   = sun_frac_total(:)   + 0.01*CTF(I_CT,I,J)*sun_frac(:)
      else if (totalCT .lt. 0) then
             !Send ERROR message!
      else   !totalCT == 0
             !Do nothing: (default values)
      endif
   
      enddo
   enddo

contains
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION Calcbeta
   !   Calculates the solar zenith angle
   !   Code originally developed by Alex Guenther in 1990s
   !   Coded into FORTRAN by Xuemei Wang
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION Calcbeta(Day, Lat, Hour)
         IMPLICIT NONE
         INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         !INTEGER :: Day
         REAL    :: Day
         REAL :: Rpi,Hour,Lat,SinDelta,CosDelta,A,B,Sinbeta,Calcbeta
         REAL,PARAMETER :: Pi = 3.14159, Rpi180 = 57.29578
   !--------------------------------------------------------------------
         SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
         CosDelta = (1 - SinDelta**2.)**0.5
   
         A = SIN(Lat / Rpi180) * SinDelta
         B = COS(Lat / Rpi180) * Cosdelta
         Sinbeta = A + B * COS(2 * Pi * (Hour - 12) / 24)
         Calcbeta = ASIN(Sinbeta) * 57.29578
   END FUNCTION Calcbeta
   ! REVISE - Delete DIstomata
   !!
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION CalcEccentricity
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION CalcEccentricity(Day)
         IMPLICIT NONE
         REAL    :: Day
         !INTEGER :: Day
         INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: CalcEccentricity
         CalcEccentricity = 1 + 0.033 * COS(2*3.14*(Day-10)/365)
   END FUNCTION CalcEccentricity
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   SUBROUTINE GaussianDist
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   subroutine GaussianDist(Distgauss, Layers)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         INTEGER,INTENT(IN) ::  Layers
         REAL,DIMENSION(Layers),INTENT(OUT) :: Distgauss
         ! local variables
         INTEGER ::  i
   !----------------------------------------------------------------
         IF (Layers .EQ. 1) THEN
           Distgauss(1)   = 0.5
         ELSEIF (Layers .EQ. 3) THEN
           Distgauss(1)   = 0.112702
           Distgauss(2)   = 0.5
           Distgauss(3)   = 0.887298
         ELSEIF (Layers .EQ. 5) THEN
           Distgauss(1)   = 0.0469101
           Distgauss(2)   = 0.2307534
           Distgauss(3)   = 0.5
           Distgauss(4)   = 0.7692465
           Distgauss(5)   = 0.9530899
         ELSE
           DO i = 1, Layers
             Distgauss(i) = (i - 0.5) / Layers
           ENDDO
         ENDIF
         RETURN
   end subroutine GaussianDist
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   SUBROUTINE SolarFractions
   !   Based on actual and potential max solar radiation:
   !   Determine the fraction of solar radiation that is 
   !   diffuse PPFD, direct PPFD, diffuse near IR, direct near IR 
   !
   !   Originally developed by Alex Guenther in 1990s
   !   Modified in 2010
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   subroutine SolarFractions( Solar, Maxsolar, Qdiffv,Qbeamv, Qdiffn, Qbeamn)
         IMPLICIT NONE
         INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL,INTENT(IN) :: Solar, Maxsolar
         REAL,INTENT(OUT) ::  Qdiffv,Qbeamv, Qdiffn, Qbeamn
         REAL :: FracDiff, PPFDfrac,PPFDdifFrac,Qv, Qn
         ! internal variables
         INTEGER :: I,J
         REAL ::  Transmis
   
         IF (Maxsolar  <= 0) THEN
           Transmis  = 0.5
         ELSEIF (Maxsolar < Solar) THEN
           Transmis  = 1.0
         ELSE
           Transmis  = Solar  / Maxsolar
         ENDIF
         !FracDiff is based on Lizaso 2005
         FracDiff    = 0.156 + 0.86/(1 + EXP(11.1*(Transmis -0.53)))
         !PPFDfrac is based on Goudrian and Laar 1994
         PPFDfrac    = 0.55 -Transmis*0.12
         !PPFDdifFrac is based on data in Jacovides 2007
         PPFDdifFrac = FracDiff *(1.06 + Transmis*0.4)
         ! Calculate  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine
         IF (PPFDdifFrac > 1.0) THEN
         PPFDdifFrac = 1.0
         ENDIF
         Qv     = PPFDfrac * Solar
         Qdiffv = Qv * PPFDdifFrac
         Qbeamv = Qv - Qdiffv
         Qn     = Solar - Qv
         Qdiffn =  Qn * FracDiff
         Qbeamn =  Qn - Qdiffn
         RETURN
   end subroutine SolarFractions
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
   !   Subroutine CanopyRad
   !
   !   Canopy light environment model
   !   Code originally developed by Alex Guenther in 1990s
   !   Coded into FORTRAN by Xuemei Wang
   !   based on Spitters et al. (1986), 
   !   Goudrian and van Laar (1994), Leuning (1997)
   !   Initial code 8-99, 
   !   modified 7-2000, 12-2001, 1-2017
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   SUBROUTINE CanopyRad(Distgauss, Layers, LAI, Sinbeta,             &
                   Qbeamv, Qdiffv, Qbeamn, Qdiffn, Cantype,      &
                   Canopychar, Sunfrac, QbAbsV, QdAbsV, QsAbsV,  &
                   QbAbsn, QdAbsn, QsAbsn, SunQv,                &
                   ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD,  &
                   NrCha, NrTyp)
   
         IMPLICIT NONE
         ! input
         INTEGER,INTENT(IN) :: Layers, NrCha, NrTyp, Cantype
         REAL,INTENT(IN) :: Qbeamv,Qdiffv,Sinbeta,LAI,Qbeamn,Qdiffn
         REAL,DIMENSION(Layers),INTENT(IN) :: Distgauss
         ! output
         REAL,INTENT(OUT) :: QbAbsV, QbAbsn
         REAL,DIMENSION(Layers),INTENT(OUT) :: ShadePPFD, SunPPFD, &
                       QdAbsv, QsAbsv, QsAbsn, ShadeQv,  SunQn,  &
                       QdAbsn, SunQv, ShadeQn, Sunfrac 
         REAL,DIMENSION(NrCha,NrTyp),INTENT(OUT) :: Canopychar
         ! internal variables
         INTEGER :: i
         REAL :: ScatV, ScatN, RefldV, RefldN, ReflbV, ReflbN,     & 
                Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, Cluster, & 
                QdAbsVL, QsAbsVL, QdAbsNL, QsAbsNL, CANTRAN, LAIadj
         
         
         REAL,PARAMETER :: ConvertShadePPFD = 4.6
         REAL,PARAMETER :: ConvertSunPPFD = 4.0
   
         ! adjust LAI for canopy transparency
         CANTRAN = Canopychar(17,Cantype)
         LAIadj = LAI / ( 1 - CANTRAN )
   
        IF (((Qbeamv  + Qdiffv ) > 0.001) .AND.  (Sinbeta  > 0.002) .AND. (LAIadj  > 0.001)) THEN       ! Daytime
   
           ! Scattering coefficients (scatV,scatN), diffuse and beam reflection 
           ! coefficients (ref..) for visible or near IR
           ScatV   = Canopychar(5,Cantype)
           ScatN   = Canopychar(6,Cantype)
           RefldV  = Canopychar(7,Cantype)
           RefldN  = Canopychar(8,Cantype)
           Cluster = Canopychar(9,Cantype)
           !        print*,'cluster',  Cluster
           ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
           Kb = Cluster * 0.5 / Sinbeta
           ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
           Kd = 0.8 * Cluster
           ! (0.8 assumes a spherical leaf angle distribution)
   
           Call CalcExtCoeff(Qbeamv,ScatV,Kb,Kd,ReflbV,KbpV,KdpV,QbAbsV)
           Call CalcExtCoeff(Qbeamn,ScatN,Kb,Kd,ReflbN,KbpN,KdpN,QbAbsn)
   
           DO i = 1,layers
             ! LAI depth at this layer
             LAIdepth   = LAIadj  * Distgauss(i)
             !fraction of leaves that are sunlit
             Sunfrac(i) = EXP(-Kb * LAIdepth)
   
             Call CalcRadComponents(Qdiffv , Qbeamv , kdpV, kbpV, kb, scatV, refldV, reflbV, LAIdepth, QdAbsVL, QsAbsVL)
             Call CalcRadComponents(Qdiffn , Qbeamn , kdpN, kbpN, kb, scatN, refldN, reflbN, LAIdepth, QdAbsNL, QsAbsNL)
   
             ShadePPFD(i) = (QdAbsVL + QsAbsVL) * ConvertShadePPFD/(1 - scatV)
             SunPPFD(i) = ShadePPFD(i) + (QbAbsV* ConvertSunPPFD/(1 - scatV))
             QdAbsV(i) = QdAbsVL
             QsAbsV(i) = QsAbsVL
             QdAbsn(i) = QdAbsNL
             QsAbsn(i) = QsAbsNL
             ShadeQv(i) = QdAbsVL + QsAbsVL
             SunQv(i)   = ShadeQv(i) + QbAbsV
             ShadeQn(i) = QdAbsNL + QsAbsNL
             SunQn(i)   = ShadeQn(i) + QbAbsn
           ENDDO
   
        ELSE                           ! Night time
          QbAbsV  = 0
          QbAbsn   = 0
          Sunfrac(:)   = 0.2
          SunQn(:)     = 0
          ShadeQn(:)   = 0
          SunQv(:)     = 0
          ShadeQv(:)   = 0
          SunPPFD(:)   = 0
          ShadePPFD(:) = 0
          QdAbsV(:)    = 0
          QsAbsV(:)    = 0
          QdAbsn(:)    = 0
          QsAbsn(:)    = 0
        ENDIF
        RETURN
   END SUBROUTINE CanopyRad
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   Subroutine CalcExtCoeff
   !   Calculate light extinction coefficients
   !   Code originally developed by Alex Guenther in 1990s
   !   Coded into FORTRAN by Xuemei Wang
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   SUBROUTINE CalcExtCoeff(Qbeam,scat,kb,kd,reflb,kbp,kdp,QbeamAbsorb)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL,INTENT(IN) :: Qbeam, scat, Kb, Kd
         REAL,INTENT(OUT) :: Reflb, Kbp, Kdp, QbeamAbsorb
         ! local variables
         REAL :: P
   
         P     = (1 - scat)**0.5
         Reflb = 1 - Exp((-2 * ((1 - P) / (1 + P)) * kb) / (1 + kb))
         ! Extinction coefficients
         Kbp   = Kb * P
         Kdp   = Kd * P
         ! Absorbed beam radiation
         QbeamAbsorb = kb * Qbeam * (1 - scat)
         RETURN
   END SUBROUTINE CalcExtCoeff
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   Subroutine CalcRadComponents
   !   Code originally developed by Alex Guenther in 1990s
   !   Coded into FORTRAN by Xuemei Wang
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb, scat, refld, reflb, LAIdepth, QdAbs, QsAbs)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL,INTENT(IN)    :: Qdiff,Qbeam,kdp,kbp,kb,scat,refld,reflb,LAIdepth
         REAL,INTENT(OUT)   :: QdAbs, QsAbs
   !-------------------------------------------------------------------
         QdAbs = Qdiff *   Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
         QsAbs = Qbeam * ((Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth)) - (Kb * (1 - Scat) * Exp(-Kb * LAIdepth)))
         RETURN
   END SUBROUTINE CalcRadComponents
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   Subroutine CanopyEB
   !
   !   Canopy energy balance model for estimating leaf temperature
   !   Coded into FORTRAN by Xuemei Wang
   !   Code developed by Alex Guenther in 1990s
   !   based on Goudrian and Laar (1994) and Leuning (1997)
   !   Initial code 8-99, modified 7-2000 and 12-2001
   !   Modified in 1-2017 by Alex Guenther and Ling Huang
   !   to correct IR balance and atmos. emissivity
   !   Note: i denotes an array containing a vertical profile 
   !         through the canopy with 0 (above canopy conditions) 
   !         plus 1 to number of canopy layers
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   
   SUBROUTINE CanopyEB(Trate, Layers, Distgauss, Canopychar,            &
                Cantype, TairK, HumidairPa, Ws,                         &
                SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn,     &
                Sunleaftk, SunleafSH, SunleafLH,                        &
                SunleafIR, Shadeleaftk, ShadeleafSH,                    &
                ShadeleafLH, ShadeleafIR, NrCha, NrTyp, Ws0,            &
                TairK0, HumidairPa0)
        IMPLICIT NONE
        INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
   
   ! inputs
        INTEGER,INTENT(IN) :: NrCha, NrTyp, Layers, Cantype
        REAL,INTENT(IN) :: Trate, TairK0, HumidairPa0, Ws0
        REAL,DIMENSION(Layers),INTENT(IN) ::  Distgauss, SunQv,ShadeQv,SunQn, ShadeQn, SunPPFD, ShadePPFD
        REAL,DIMENSION(NrCha, NrTyp),INTENT(IN)  :: Canopychar
   
   ! outputs
         REAL,DIMENSION(Layers),INTENT(OUT) :: HumidairPa, Ws, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, TairK, Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR
   
   ! local variables
         INTEGER :: i
         REAL :: Cdepth, Lwidth, Llength, Cheight, Eps, TranspireType, Deltah, EmissAtm, LeafIR, IRin,IRout
   !     &         Deltah, UnexposedLeafIRin, ExposedLeafIRin, IRin,IRout
         REAL,DIMENSION(Layers) :: Ldepth, Wsh
   !-----------------------------------------------------------------------
   
         Cdepth        = Canopychar(1, Cantype)
         Lwidth        = Canopychar(2, Cantype)
         Llength       = Canopychar(3, Cantype)
         Cheight       = Canopychar(4, Cantype)
         Eps           = Canopychar(10,Cantype)
         TranspireType = Canopychar(11,Cantype)
   
         IF (TairK0  > 288) THEN
   ! Pa m-1  (humidity profile for T < 288)
           Deltah =  Canopychar(14, Cantype) / Cheight
         ELSEIF (TairK0  > 278) THEN
           Deltah =(Canopychar(14,Cantype)-((288-TairK0)/10) * (Canopychar(14,Cantype)-Canopychar(15,Cantype)))/Cheight
         ELSE
   ! Pa m-1  (humidity profile for T <278)
           Deltah = Canopychar(15, Cantype) / Cheight
         ENDIF
   
         Ldepth(:)     = Cdepth * Distgauss(:)
         TairK(:)      = TairK0  + (Trate  * Ldepth(:))      ! check this
         HumidairPa(:) = HumidairPa0  + (Deltah * Ldepth(:))
   
         Wsh(:) = (Cheight-Ldepth(:)) - (Canopychar(16,Cantype) * Cheight)
         Ws(:)  = (Ws0*LOG(Wsh(:))/LOG(Cheight-Canopychar(16,Cantype) * Cheight))
         WHERE (Wsh(:) < 0.001) Ws(:) = 0.05
   
         DO i=1,Layers
   
            ! REVISE - Replace UnexposedLeafIR with LeafIR
            
            !        IRin     = UnexposedLeafIRin(TairK(i), Eps)
            !        ShadeleafIR(i) = 2 * IRin
            !        SunleafIR(i) = 0.5*ExposedLeafIRin(HumidairPa0,TairK0)+1.5*IRin
            
            ! Apparent atmospheric emissivity for clear skies: 
            ! function of water vapor pressure (Pa) 
            ! and ambient Temperature (K) based on Brutsaert(1975) 
            ! referenced in Leuning (1997)
            EmissAtm        = 0.642 * (HumidairPa(i) / TairK(i))**(1./7.)   
            IRin            = LeafIR (TairK(i), EmissAtm)
            ShadeleafIR(i)  = IRin
            SunleafIR(i)    = IRin
   
         ! Sun
           CALL LeafEB(SunPPFD(i), SunQv(i) + SunQn(i),                    &
                      SunleafIR(i), Eps, TranspireType, Lwidth, Llength,   &
                      TairK(i), HumidairPa(i), Ws(i),                      &
                      Sunleaftk(i), SunleafSH(i),SunleafLH(i),             &
                      IRout )
   
            SunleafIR(i) = SunleafIR(i) - IRout
   
         ! Shade
           CALL LeafEB(ShadePPFD(i), ShadeQv(i)+ShadeQn(i),                &
                        ShadeleafIR(i),Eps,TranspireType, Lwidth,Llength,  &
                        TairK(i), HumidairPa(i), Ws(i),                    &
                        Shadeleaftk(i), ShadeleafSH(i),ShadeleafLH(i),     &
                        IRout)
   
            ShadeleafIR(i) = ShadeleafIR(i) - IRout
         ENDDO
   
         RETURN
   END SUBROUTINE CanopyEB
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   Subroutine LeafEB
   !
   !   Leaf energy balance
   !   Code originally developed by Alex Guenther in 1990s
   !   Coded into FORTRAN by Xuemei Wang
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType,         &
            Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf,      &
            SH, LH, IRout)
   
         IMPLICIT NONE
         INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
   
         REAL,INTENT(IN) :: Eps, TranspireType, Lwidth, Llength,PPFD, Q, IRin, TairK, HumidairPa, Ws
         REAL,INTENT(OUT) :: IRout, Tleaf, SH, LH
   
   ! local variables
   
         INTEGER :: i
         REAL :: HumidAirKgm3,GHforced,StomRes,IRoutairT,LHV,LatHv,Ws1,      &
           LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIR,                 &
           GH1,SH1,LH1,E1,ResSC,IRout1,GH !,ConvertHumidityPa2kgm3
   !     &        LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIRout,   
   !----------------------------------------------------
   
         IF (Ws <= 0) THEN
           Ws1 = 0.001
         ELSE
           Ws1 = Ws
         ENDIF
   
         ! Air vapor density kg m-3
         HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)
   
         ! Heat convection coefficient (W m-2 K-1) for forced convection. 
         ! Nobel page 366
         GHforced = 0.0259 / (0.004 * ((Llength / Ws)**0.5))
   
         ! Stomatal resistence s m-1
         StomRes  = ResSC(PPFD)
   
         ! REVISE - Replace LeafIRout with LeafIR
         !      IRoutairT = LeafIROut(tairK, eps)
         !XJ      IRoutairT  = LeafIR(TairK + Tdelt, Eps)
          IRoutairT = LeafIR(TairK, Eps)
   
         ! Latent heat of vaporization (J Kg-1)
         LatHv = LHV(TairK)
   
         ! Latent heat flux
         LHairT = LeafLE(TairK,HumidAirKgm3,LatHv,GHforced,StomRes, TranspireType)
   
         E1 = (Q + IRin - IRoutairT - LHairT)
         IF (E1 .EQ. 0.) THEN
           E1 = -1.
         ENDIF
   
         Tdelt = 1.
         Balance = 10.
         DO i = 1, 10
           IF (ABS(Balance) > 2) THEN
             ! Boundary layer conductance
             GH1 = LeafBLC(GHforced, Tdelt, Llength)
             ! Convective heat flux
             SH1 = LeafH(Tdelt, GH1)
             ! Latent heat of vaporization (J Kg-1)
             LatHv = LHV(TairK + Tdelt)
             LH = LeafLE(TairK + Tdelt, HumidAirKgm3, LatHv, GH1, StomRes, TranspireType)
             LH1 = LH - LHairT
             ! REVISE - Replace LeafIROut with LeafIR
             !          IRout  = LeafIROut(TairK + Tdelt, Eps)
             IRout  = LeafIR(TairK + Tdelt, Eps)
             IRout1 = IRout - IRoutairT
             Tdelt  = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
             Balance = Q + IRin - IRout - SH1 - LH
           ENDIF
         ENDDO
   
         If (Tdelt > 10)  Tdelt = 10
         If (Tdelt < -10) Tdelt = -10
   
         Tleaf = TairK + Tdelt
         GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
         SH    = LeafH(Tleaf - TairK, GH)
         LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, GH, StomRes, TranspireType)
   
         ! REVISE - Replace LeafIROut with LeafIR
         !      IRout = LeafIROut(Tleaf, Eps)
         IRout = LeafIR(Tleaf, Eps)
   
         RETURN
   END SUBROUTINE LeafEB
   
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION WaterVapPres
   !
   !   Convert water mixing ratio (kg/kg) to water vapor pressure 
   !   (Pa or Kpa depending on units of input )
   !   Mixing ratio (kg/kg), temp (C), pressure (KPa)
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION WaterVapPres(Dens, Pres, WaterAirRatio)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Dens, Pres, WaterVapPres, WaterAirRatio
         WaterVapPres = (Dens / (Dens + WaterAirRatio)) * Pres
   END FUNCTION WaterVapPres
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION Stability
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION Stability(Canopychar, Cantype, Solar, NrCha, NrTyp)
         IMPLICIT NONE
         INTEGER :: Cantype, NrCha, NrTyp
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Solar, Trateboundary, Stability
         REAL,DIMENSION(NrCha, NrTyp)  :: Canopychar
   
         Trateboundary = 500
         IF (Solar > Trateboundary) THEN
               ! Daytime temperature lapse rate
           Stability = Canopychar(12, Cantype)
         ELSEIF (Solar > 0) THEN
           Stability = Canopychar(12, Cantype) - ((Trateboundary - Solar) / Trateboundary) * (Canopychar(12, Cantype) - Canopychar(13, Cantype))
         ELSE
               ! Nightime temperature lapse rate
            Stability = Canopychar(13, Cantype)
         ENDIF
   END FUNCTION Stability
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION ConvertHumidityPa2kgm3
   !
   !   Saturation vapor density  (kg/m3)
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION ConvertHumidityPa2kgm3(Pa, Tk)
         implicit none
         integer,parameter :: real_x = selected_real_kind(p=14, r=30)
         REAL              :: ConvertHumidityPa2kgm3, Pa, Tk
         ConvertHumidityPa2kgm3 = 0.002165 * Pa / Tk
   END FUNCTION ConvertHumidityPa2kgm3
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION ResSC
   !
   !   Leaf stomatal cond. resistance s m-1
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION ResSC(Par)
         IMPLICIT NONE
         INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Par, SCadj, ResSC
         SCadj = ((0.0027 * 1.066 * Par) / ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))
         IF (SCadj < 0.1) THEN
           ResSC = 2000
         ELSE
           ResSC = 200 / SCadj
         ENDIF
   END FUNCTION ResSC
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION LeafIR
   !
   !   Calculate IR transfer between leaf and air
   !   Added by Alex Guenther and Ling Huang to replace previous
   !   MEGAN2.1 IR balance functions
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION LeafIR(Tk, Eps)
          IMPLICIT NONE
          INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
          REAL :: Eps, Tk, LeafIR
          LeafIR = Eps * Sb * (2 * (Tk**4.))
   END FUNCTION LeafIR
   
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION LHV
   !
   !   Latent Heat of vaporization(J Kg-1) from Stull p641
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION LHV(Tk)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Tk, LHV
         ! REVISE - Replace 273 with 273.15
         !      LHV = 2501000 - (2370 * (Tk - 273))
         LHV = 2501000 - (2370 * (Tk - 273.15))
   END FUNCTION LHV
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION LeafLE
   !
   !   Latent energy term in Energy balance
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType, SvdTk,LeafRes, Vapdeficit, LeafLE, LE
         LeafRes    = (1 / (1.075 * (GH / 1231))) + StomRes
         Vapdeficit = (SvdTk(Tleaf) - Ambvap)
         ! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / 
         !                 leaf resistence (s m-1)
         LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit
         IF (LE < 0) THEN
           LeafLE = 0
         ELSE
           LeafLE = LE
         ENDIF
   END FUNCTION  LeafLE
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION LeafBLC
   !
   !   Boundary layer conductance
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION LeafBLC(GHforced, Tdelta, Llength)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: GHforced, Tdelta, Llength, Ghfree, LeafBLC
         !----------------------------------------------------------------
         ! This is based on Leuning 1995 p.1198 except using molecular 
         ! conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
         ! diffusivity so that you end up with a heat convection coefficient 
         ! (W m-2 K-1) instead of a conductance for free convection
         IF (Tdelta >= 0) THEN
            GhFree = 0.5 * 0.00253 * ((160000000 * Tdelta / (Llength**3.))**0.25) / Llength
         ELSE
           GhFree = 0
         ENDIF
         LeafBLC = GHforced + GhFree
   END FUNCTION LeafBLC
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION LeafH
   !
   !   Convective energy term in Energy balance (W m-2 heat flux 
   !      from both sides of leaf)
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   FUNCTION LeafH(Tdelta, GH)
         IMPLICIT NONE
         INTEGER,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Tdelta, GH, LeafH
         LeafH = 2 * GH * Tdelta! 2 sides X conductance X Temperature gradient
   END FUNCTION LeafH
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   !
   !   FUNCTION SvdTk
   !
   !   Saturation vapor density  (kg/m3)
   !
   !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   Function SvdTk(Tk)
         IMPLICIT NONE
         INTEGER ,PARAMETER :: real_x = SELECTED_REAL_KIND(p=14, r=30)
         REAL :: Tk, Svp, SvdTk
         ! Saturation vapor pressure (millibars)
         Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)
         SvdTk = 0.2165 * Svp / Tk
   END FUNCTION  SvdTk
end subroutine megcan

end module megcanopy
