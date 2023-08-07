module SCRIP
!  This is a module that tries to encapsulate SCRIP procedures

   use SCRIP_KindsMod             ! module defining data types
   use SCRIP_CommMod              ! for initializing comm environment
   use SCRIP_ErrorMod             ! SCRIP error checking and logging
   use SCRIP_IOUnitsMod           ! manages I/O units
   use SCRIP_ConfigMod            ! SCRIP configuration module
   use SCRIP_InitMod              ! SCRIP initialization
   use constants                  ! module for common constants
   use timers                     ! CPU timers
   use grids                      ! module with grid information
   use remap_vars                 ! common remapping variables
   use remap_conservative         ! routines for conservative remap
   use remap_distance_weight      ! routines for dist-weight remap
   use remap_bilinear             ! routines for bilinear interp
   use remap_bicubic              ! routines for bicubic  interp
   use remap_write                ! routines for remap output

contains


subroutine remap(g,iFile,oFile,varlist, mapMethod)

   implicit none 
   !-----------------------------------------------------------------------
   ! external variables
   !-----------------------------------------------------------------------
   type(grid_type) :: g                !grid/proj specs of dest grid (will be saved on oFile)
   character(254)  :: iFile, oFile     !input and output File.
   character(50)   :: varlist(:)       !list of variable names to remap on iFile.
   integer         :: mapMethod        !bilinear, bicubic, conserv (1st and 2nd).
   !-----------------------------------------------------------------------
   !  local variables
   !-----------------------------------------------------------------------
   integer (SCRIP_i4) :: errorCode      ! error flag
   
   character (SCRIP_charLength) :: &
      gridFile1,    &! filename of grid file containing grid1
      gridFile2,    &! filename of grid file containing grid2
      interpFile1,  &! filename for output remap data (map1)
      interpFile2,  &! filename for output remap data (map2)
      mapName1,     &! name for mapping from grid1 to grid2
      mapName2,     &! name for mapping from grid2 to grid1
      mapMethod,    &! choice for mapping method
      normalizeOpt, &! option for normalizing weights
      outputFormat   ! option for output conventions
   
   character (12), parameter :: rtnName = 'SCRIP_driver'
   
   integer (SCRIP_i4) :: n, iunit  ! dummy counter & unit number for input configuration file
   
   !-----------------------------------------------------------------------
   !  initialize communication environment and SCRIP package
   !-----------------------------------------------------------------------
   errorCode = SCRIP_Success
   
   call SCRIP_CommInitMessageEnvironment
   call SCRIP_Initialize(errorCode)
   !if (SCRIP_errorCheck(errorCode, rtnName, 'error initializing SCRIP')) &
   !   call SCRIP_driverExit(errorCode)
   
   !-----------------------------------------------------------------------
   !  initialize timers
   !-----------------------------------------------------------------------
      call timers_init
      do n=1,max_timers
         call timer_clear(n)
      end do
   !-----------------------------------------------------------------------
   !  read input namelist
   !-----------------------------------------------------------------------
                                       !&remapInputs
    num_maps=1                         ! num_maps = 2                      !number of mappings to be computed.(1: forwards, 2: both directions)
    mapMethod=mapMethod                ! mapMethod = ’conservative’        !method
                                                                                                                                                      
    gridFile1=iFile                    ! gridFile1 = ’grid_1_file_name’         !
    gridFile2=oFile                    ! gridFile2 = ’grid_2_file_name’         !
    interpFile1="interp_"//trim(iFile) ! interpFile1 = ’map_1_output_file_name’ !
    interpFile2="interp_"//trim(oFile) ! interpFile2 = ’map_2_output_file_name’ ! 
    mapName1='mapName1'                ! mapName1 = ’name_for_mapping_1’        !
    mapName2='mapName2'                ! mapName2 = ’name_for_mapping_2’        !
                                       
    outputFormat='scrip'               ! outputOpt = ’scrip’                    !format of output file. 
                                       
    restrict_type='latitude'           ! restrict_type = ’latitude’        !'latitude' or 'latitude-longitude' search restriction to avoid N^2 search
    num_srch_bins=90                   ! num_srch_bins = 90                !number of bins for search
    normalizeOpt='fracArea'            ! normalizeOpt = ’fracArea’         !'fracArea', 'none', 'dstArea' (normalization for conservative remaping)
    luse_grid1_area=.false.            ! luse_grid1_area = .false.
    luse_grid2_area=.false.            ! luse_grid2_area = .false.
    npseg=11                           ! num_polar_segs = 11
    north_thresh=1.5_SCRIP_r8          ! npole_threshold = 1.5
    south_thresh=-2.0_SCRIP_r8         ! spole_threshold = -1.5
    nthreads=2                         ! num_threads = 2
                                       !/

    select case(mapMethod)
    case ('conservative')
       map_type = map_type_conserv
       luse_grid_centers = .false.
    case ('bilinear')
       map_type = map_type_bilinear
       luse_grid_centers = .true.
    case ('bicubic')
       map_type = map_type_bicubic
       luse_grid_centers = .true.
    case ('distwgt')
       map_type = map_type_distwgt
       luse_grid_centers = .true.
    case default
       !call SCRIP_ErrorSet(errorCode, rtnName, 'unknown mapping method')
       !call SCRIP_driverExit(errorCode)
    end select

    select case(trim(normalizeOpt))
    case ('none')
       norm_opt = norm_opt_none
    case ('fracArea')
       norm_opt = norm_opt_frcarea
    case ('destArea')
       norm_opt = norm_opt_dstarea
    case default
       !call SCRIP_ErrorSet(errorCode, rtnName, 'unknown normalization option')
       !call SCRIP_driverExit(errorCode)
    end select
    
!-----------------------------------------------------------------------
!  initialize grid information for both grids
!-----------------------------------------------------------------------

   call grid_init(gridFile1, gridFile2, errorCode)

   !if (SCRIP_ErrorCheck(errorCode, rtnName, 'Error initializing grids')) &
   !   call SCRIP_driverExit(errorCode)

   write(SCRIP_stdout, *) 'Computing remappings between: ',grid1_name
   write(SCRIP_stdout, *) '                         and  ',grid2_name

!-----------------------------------------------------------------------
!  initialize some remapping variables.
!-----------------------------------------------------------------------
   call init_remap_vars

!-----------------------------------------------------------------------
!  call appropriate interpolation setup routine based on type of
!  remapping requested.
!-----------------------------------------------------------------------
   select case(map_type)
   case(map_type_conserv)
      call remap_conserv
   case(map_type_bilinear)
      call remap_bilin
   case(map_type_distwgt)
      call remap_distwgt
   case(map_type_bicubic)
      call remap_bicub
   case default
   !   call SCRIP_ErrorSet(errorCode, rtnName, 'Invalid Map Type')
   !   call SCRIP_driverExit(errorCode)
   end select                                     

!-----------------------------------------------------------------------
!  reduce size of remapping arrays and write remapping info to file.
!-----------------------------------------------------------------------

   if (num_links_map1 /= max_links_map1) then
        call resize_remap_vars(1, num_links_map1-max_links_map1)
   endif
   if ((num_maps > 1) .and. (num_links_map2 /= max_links_map2)) then
        call resize_remap_vars(2, num_links_map2-max_links_map2)
   endif

   call write_remap( mapName1   , mapName2   , &
                     interpFile1, interpFile2, &
                    outputFormat, errorCode    )

!-----------------------------------------------------------------------
!  All done, exit gracefully
!-----------------------------------------------------------------------
   call SCRIP_driverExit(errorCode)


end subroutine

end module remap_SCRIP
