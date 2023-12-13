module mgn2mech
use megan_ini
    implicit none
     REAL, PARAMETER :: nmol2mol   = 1E-9    ! convert nanomoles to moles
contains
  subroutine convert2mech(ncols,nrows,  &
            efmaps,inper,outer)
     implicit none
     !input vars:
     integer,intent(in) :: ncols, nrows
     real   ,intent(in) :: efmaps(ncols,nrows,nclass)
     real   ,intent(in) :: inper(ncols,nrows,nclass)
     !output vars:
     real   ,intent(inout) :: outer(ncols,nrows,nmgnspc)
     ! local variables and their descriptions:
     INTEGER :: ios,i_NO
     INTEGER :: t, s, I, N                   ! Counters
     INTEGER :: nmpmg, nmpsp, nmpmc          ! Counters
  
     REAL, ALLOCATABLE :: tmper(:,:,:)       ! Temp emission buffer
     !REAL :: GAMNO(ncols,nrows)         ! NO emission factor
     !REAL :: BDSNP_NO(ncols,nrows)      ! NO emissions for BDSNP algorithm (nanomol/m^2/s)
     !REAL, INTENT(IN)     ::  inper       (19, NCOLS, NROWS )
     !REAL, INTENT(IN)     ::  NO_IN       (NCOLS, NROWS )
     !REAL, INTENT(OUT)    ::  outer       (NMGNSPC, NCOLS, NROWS ) ! CB6
  
     !.....2) Conversion from MGN 20 to speciated 201
     !-----------------------------------------------------------------------
     !...  Allocate memory
     ALLOCATE ( tmper( ncols, nrows, n_spca_spc), STAT = ios )
  
     i_NO = 8 ! this was 20 for megan 3.1 
  
     tmper = 0.
     outer = 0.
  
     !@IF ( .NOT. BDSNP_MEGAN ) THEN
     !@  GAMNO = NO_IN
     !@ELSE
     !@  BDSNP_NO = NO_IN
     !@ENDIF
  
     DO s = 1, N_SMAP_SPC
       nmpmg = mg20_map(s)
       nmpsp = spca_map(s)
       !@IF ( nmpmg .NE. i_NO ) then !...  Not NO
         tmper(:,:,nmpsp) = inper(:,:,nmpmg) * efmaps(:,:,nmpmg)  * effs_all(s)
       !@ELSEIF ( nmpmg .EQ. i_NO ) then
         !@!!-----------------NO Stuff-----------------------
         !@IF ( .NOT. BDSNP_MEGAN ) THEN
         !@!     GAMNO is emission activity factor
         !@   tmper(:,:,nmpsp) = GAMNO(:,:) * efmaps(:,:,i_NO)        * effs_all(s)
         !@ELSE
         !@! directly use BDSNP soil NO
         !@  tmper(nmpsp,:,:) = BDSNP_NO(:,:)
         !@ENDIF
         !@!-----------------end of NO----------------------
       !@ENDIF     !IF ( nmpmg .NE. i_NO ) then
     ENDDO ! End species loop
     !-----------------------------------------------------------------------
     !.....3) Conversion from speciated species to MECHANISM species
     !-----------------------------------------------------------------------
      DO s = 1, n_spca_spc
         tmper(s,:,:) = tmper(s,:,:) * nmol2mol
      ENDDO
  
      ! lumping to MECHANISM species
      DO s = 1, n_scon_spc
        nmpsp = spmh_map(s)         ! Mapping value for SPCA
        nmpmc = mech_map(s)         ! Mapping value for MECHANISM
  
        IF ( nmpmc .NE. 999 ) THEN
           outer(:,:,nmpmc) = outer(:,:,nmpmc) +  (tmper(:,:,nmpsp) * conv_fac(s))
        ENDIF
      ENDDO ! End species loop
  
  end subroutine convert2mech
end module mgn2mech
