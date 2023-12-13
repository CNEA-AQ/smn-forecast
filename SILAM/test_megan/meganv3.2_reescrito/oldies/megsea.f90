module meg_sea
      INTEGER, PARAMETER :: NrTyp = 6          ! Number of canopy types
      REAL,PARAMETER :: d1 = 0.04

      !-- WWLT is wilting point (M^3/M^3) (JN90)
      REAL WWLT_JN90(16)
      DATA WWLT_JN90 /  0.068, 0.075, 0.114, 0.179, 0.155, 0.175, 0.218, 0.250, 0.219, 0.283, 0.286, 0.286, 0.286, 0.286, 0.286, 0.286 /

     !-- WWLT is wilting point (M^3/M^3)
     REAL  WWLT_NOAH(19)
     DATA  WWLT_NOAH / 0.024, 0.057, 0.081, 0.123, 0.064, 0.128, 0.168, 0.212, 0.196, 0.239, 0.264, 0.285, 0.118, 0.066, 0.009, 0.049, 0.264, 0.009, 0.015 /

contains

subroutine megsea(              &
             NCOLS,NROWS,       &
             LSM, SLTYP, SOILM, &
             GAMMA_SM           )

      implicit none
     ! input  variables
     integer, intent(in)           :: ncols,nrows
     character(len=4),intent(in)   :: LSM          !land surface model 
     integer, intent(in)           :: sltyp(NCOLS,NROWS)     !soil type     [lsm_id]
     real,    intent(in)           :: soilm(NCOLS,NROWS)     !soil moisture [1]
     ! output variables
     real, intent(out)             :: gamma_sm(ncols, nrows) ! Soil moisture activity for isoprene
    
     ! local vars:
     integer :: i,j 
     real, allocatable :: wwlt(:)
     real :: t1, wilt
      SELECT CASE (LSM)
             CASE ('NOAH')
                allocate(wwlt(size(wwlt_noah)));wwlt=wwlt_noah;
             CASE ('JN90')
                allocate(wwlt(size(wwlt_jn90)));wwlt=wwlt_jn90;
             CASE DEFAULT
                allocate(wwlt(size(wwlt_jn90)));wwlt=wwlt_jn90;
      END SELECT

     do j = 1, nrows
        do i = 1, ncols

         wilt = wwlt(sltyp(i,j))
         t1 = wilt + d1
         IF ( soilm(i,j) < wilt ) THEN
             gamma_sm(i,j) = 0

         else if ( soilm(i,j) >= wilt .and. soilm(i,j) < t1 ) then
             gamma_sm(i,j) = (soilm(i,j) - wilt)/d1

         else
             gamma_sm(i,j) = 1
         end if

       end do 
     end do 
       
end subroutine megsea

end module meg_sea
