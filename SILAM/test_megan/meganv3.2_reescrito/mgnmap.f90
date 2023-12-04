module megan_ini
   !
   ! This module determines the speciation map files to use and populates 
   ! the MEGAN_NAMES array. This array is used in EMIS_DEFN.F via the 
   ! variable vdemis_me in the MIOG stream.
   !
   implicit none
   !Number of emission classes
   INTEGER, PARAMETER :: NCLASS = 19
   INTEGER, PARAMETER :: NEMIS  = NCLASS
   ! number of emission classes

    integer, save :: nmgnspc           !number of megan     species
    integer, save :: n_scon_spc        !number of mechanism species
    character( 16 ), allocatable :: megan_names(:)     ! megan species names
    integer, allocatable ::  spmh_map(:),mech_map(:)   ! speciated species name

    real, allocatable :: conv_fac(:)        
    real,allocatable :: mech_mwt(:)
    character( 16 ), allocatable :: mech_spc(:)

    INCLUDE 'tables/SPC_NOCONVER.EXT'
    INCLUDE 'tables/SPC_CB05.EXT'
    INCLUDE 'tables/SPC_CB6.EXT'
    INCLUDE 'tables/SPC_CB6_AE7.EXT'
    INCLUDE 'tables/SPC_RACM2.EXT'        ! new in MEGAN3
    INCLUDE 'tables/SPC_CRACMM.EXT'       ! new in CMAQ 5.4
    INCLUDE 'tables/MAP_CV2CB05.EXT'
    INCLUDE 'tables/SPC_SAPRC07.EXT'      ! new in MEGAN3
    INCLUDE 'tables/SPC_SAPRC07T.EXT'     ! new in MEGAN3
    INCLUDE 'tables/MAP_CV2CB6.EXT'      
    INCLUDE 'tables/MAP_CV2CB6_AE7.EXT'      
    INCLUDE 'tables/MAP_CV2RACM2.EXT'
    INCLUDE 'tables/MAP_CV2CRACMM.EXT'
    INCLUDE 'tables/MAP_CV2SAPRC07.EXT'
    INCLUDE 'tables/MAP_CV2SAPRC07T.EXT'

   contains

   subroutine megan_map(mechanism)

     implicit none
       character( 16 )    :: mechanism     ! mechanism name
       integer, parameter :: nmechs = 16   ! dimension for number of mechanisms considered
       integer ios,indx
       integer i

       TYPE MIOG_MECH_TYPE
            CHARACTER( 32 ) :: CHEMMECH
            CHARACTER( 16 ) :: MIOGMECH
       END TYPE MIOG_MECH_TYPE

       TYPE(MIOG_MECH_TYPE) :: MIOG_MECH_MAP( NMECHS ) = [       &
          miog_mech_type('CB05E51_AE6_AQ         ','CB05    '),  &
          miog_mech_type('CB05EH51_AE6_AQ        ','CB05    '),  &
          miog_mech_type('CB05MP51_AE6_AQ        ','CB05    '),  &
          miog_mech_type('CB05TUCL51_AE6_AQ      ','CB05    '),  &
          miog_mech_type('CB6R3_AE6_AQ           ','CB6     '),  &
          miog_mech_type('CB6MP_AE6_AQ           ','CB6     '),  &
          miog_mech_type('CB6R5HAP_AE7_AQ        ','CB6_ae7 '),  &
          miog_mech_type('CB6R3_AE7_AQ           ','CB6_ae7 '),  &
          miog_mech_type('CB6R5_AE7_AQ           ','CB6_ae7 '),  &
          miog_mech_type('CB6R5M_AE7_AQ          ','CB6_ae7 '),  &
          miog_mech_type('RACM2_AE6_AQ           ','RACM2   '),  &
          miog_mech_type('SAPRC07TC_AE6_AQ       ','SAPRC07T'),  &
          miog_mech_type('SAPRC07TIC_AE7I_AQ     ','SAPRC07T'),  &
          miog_mech_type('SAPRC07TIC_AE7I_AQKMT2 ','SAPRC07T'),  &   
          miog_mech_type('CRACMM1_AQ             ','CRACMM  '),  &
          miog_mech_type('CRACMM1AMORE_AQ        ','CRACMM  ')   ]

          !INDX = INDEX1( MECHNAME, NMECHS, MIOG_MECH_MAP%CHEMMECH )
          !MECHANISM = MIOG_MECH_MAP( INDX )%MIOGMECH
 
          SELECT CASE ( TRIM(MECHANISM) )
            CASE ('SAPRC07')
              n_scon_spc = n_saprc07
              NMGNSPC = n_saprc07_spc
            CASE ('SAPRC07T')
              n_scon_spc = n_saprc07t
              NMGNSPC = n_saprc07t_spc
            CASE ('CB05')
              n_scon_spc = n_cb05
              NMGNSPC = n_cb05_spc
            CASE ('CB6')
              n_scon_spc = n_cb6  ! 145
              NMGNSPC = n_cb6_spc ! 34
            CASE ('RACM2')
              n_scon_spc = n_racm2
              NMGNSPC = n_racm2_spc
            CASE ('CB6_ae7')
              n_scon_spc = n_cb6_ae7
              NMGNSPC = n_cb6_ae7_spc
            CASE ('CRACMM')
              n_scon_spc = n_cracmm
              NMGNSPC = n_cracmm_spc
            CASE DEFAULT
              print*,"Mechanism," // TRIM( MECHANISM) // ", is not identified."
          ENDSELECT
   
          allocate(MEGAN_NAMES(NMGNSPC), stat=ios)
          allocate(spmh_map(n_scon_spc), stat=ios)
          allocate(mech_map(n_scon_spc), stat=ios)
          allocate(conv_fac(n_scon_spc), stat=ios)
          allocate(mech_spc(NMGNSPC )  , stat=ios)
          allocate(mech_mwt(NMGNSPC )  , stat=ios)
  
          SELECT CASE ( TRIM(MECHANISM) )
  
            CASE ('CB05')
              spmh_map(1:n_scon_spc) = spmh_map_cb05(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_cb05(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_cb05(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_cb05(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_cb05(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =      mech_spc(1:NMGNSPC)
            CASE ('CB6')
              spmh_map(1:n_scon_spc) = spmh_map_cb6(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_cb6(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_cb6(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_cb6(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_cb6(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =     mech_spc(1:NMGNSPC)
            CASE ('RACM2')
              spmh_map(1:n_scon_spc) = spmh_map_racm2(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_racm2(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_racm2(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_racm2(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_racm2(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =       mech_spc(1:NMGNSPC)
            CASE ('SAPRC07')
              spmh_map(1:n_scon_spc) = spmh_map_saprc07(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_saprc07(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_saprc07(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_saprc07(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_saprc07(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =         mech_spc(1:NMGNSPC)
            CASE ('SAPRC07T')
              spmh_map(1:n_scon_spc) = spmh_map_saprc07t(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_saprc07t(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_saprc07t(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_saprc07t(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_saprc07t(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =          mech_spc(1:NMGNSPC)
            CASE ('CB6_ae7')
              spmh_map(1:n_scon_spc) = spmh_map_cb6_ae7(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_cb6_ae7(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_cb6_ae7(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_cb6_ae7(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_cb6_ae7(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =         mech_spc(1:NMGNSPC)
            CASE ('CRACMM')
              spmh_map(1:n_scon_spc) = spmh_map_cracmm(1:n_scon_spc)
              mech_map(1:n_scon_spc) = mech_map_cracmm(1:n_scon_spc)
              conv_fac(1:n_scon_spc) = conv_fac_cracmm(1:n_scon_spc)
              mech_spc(1:NMGNSPC)    = mech_spc_cracmm(1:NMGNSPC)
              mech_mwt(1:NMGNSPC)    = mech_mwt_cracmm(1:NMGNSPC)
              MEGAN_NAMES(1:NMGNSPC) =        mech_spc(1:NMGNSPC)
            CASE DEFAULT 
              print*,"Mapping for Mechanism," // TRIM( MECHANISM)// ", is unspecified."
          ENDSELECT
   end subroutine megan_map
end module megan_ini
