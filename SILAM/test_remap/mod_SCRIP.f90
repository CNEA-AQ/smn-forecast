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
  use remap_read                 ! routines for reading remap files
  use remap_mod                  ! module containing remapping routines
      

   use netcdf 
   use proj_silam               !proj_trans
   !private

   !public grid2gridFile           
   !public createInterpFile
   !public applyRemap
   !public gridFile2grid2d

   implicit none
   !real, parameter :: pi=3.141593
   !real, parameter :: deg2rad=pi/180.0
   type regular_grid_type
     character(12)  :: gridName             !grid-name
     character(200) :: proj4                !any proj srs string descriptor
     integer        :: nx,ny,nz             !number of cells in x-y direction (ncols, nrows, nlevs)
     real           :: dx,dy                !x-y cell dimension (x_cell, y_cell)
     real           :: xmin,ymin,xmax,ymax  !bbox and center coorindates (in the proj4 system)
   end type regular_grid_type


contains
subroutine check(status)
  integer, intent(in) :: status
  if (status /= nf90_noerr) then
    write(*,*) nf90_strerror(status)
    stop 'netcdf error'
  end if
end subroutine check

subroutine grid2gridFile(g1,g2)
  !this subroutine creates a gridFile (necesary to pass to SCRIP to build a interp_file between two grids).
  !it gets two "grid_types" with them grid and proj parameters, and build the gridFile.
  implicit none
  type(regular_grid_type) :: g1   !source grid (regular on arbitrary coordinate system)
  type(regular_grid_type) :: g2   !destination grid (silam latlon grid)

  character(100) :: gridFile
  real(8), allocatable :: center_x(:), center_y(:)
  real(8), allocatable :: corner_x(:,:), corner_y(:,:)
  integer, allocatable :: imask(:)
  integer :: grid_size, grid_rank,grid_corners
  integer :: i,j,idx
  integer :: ncid, size_dim_id,rank_dim_id,corn_dim_id,var_id
  logical :: file_exist

  gridFile = "grd_"//trim(g1%gridName)//".nc"

  !First check if file already was created!
  inquire(file=trim(gridFile),exist=file_exist)
  if( file_exist ) then; 
     print*,"File: "//trim(gridFile)//" already exist!!";
     continue;
  else
     grid_size=g1%nx*g1%ny
     grid_rank=2
     grid_corners=4

     allocate(center_x(grid_size))
     allocate(center_y(grid_size))
     allocate(corner_x(grid_corners,grid_size))
     allocate(corner_y(grid_corners,grid_size))
     allocate(imask(grid_size))

     !compute cells center
     idx=1
     do j=0,g1%ny-1
       do i=0,g1%nx-1
          center_x(idx)=g1%xmin+0.5*g1%dx+g1%dx*i
          center_y(idx)=g1%ymin+0.5*g1%dy+g1%dy*j
          idx=idx+1
       enddo
     enddo
     !compute cells corner coordinates:
     corner_x(1,:)=center_x(:)-g1%dx       
     corner_x(2,:)=center_x(:)-g1%dx       ! 2---------4
     corner_x(3,:)=center_x(:)+g1%dx       ! |         |
     corner_x(4,:)=center_x(:)+g1%dx       ! |    o    |
     corner_y(1,:)=center_y(:)-g1%dy       ! |         |
     corner_y(2,:)=center_y(:)+g1%dy       ! 1---------3
     corner_y(3,:)=center_y(:)-g1%dy      
     corner_y(4,:)=center_y(:)+g1%dy                     
     !transform to latlon
     if (g1%proj4 /= g2%proj4) then
        call proj_trans(g1%proj4, g2%proj4, center_x     , center_y     , grid_size)     !  2---------------4 
        call proj_trans(g1%proj4, g2%proj4, corner_x(1,:), corner_y(1,:), grid_size)     !   \             /
        call proj_trans(g1%proj4, g2%proj4, corner_x(2,:), corner_y(2,:), grid_size)     !    \     o     /
        call proj_trans(g1%proj4, g2%proj4, corner_x(3,:), corner_y(3,:), grid_size)     !     \         /
        call proj_trans(g1%proj4, g2%proj4, corner_x(4,:), corner_y(4,:), grid_size)     !      1-------3
     endif
     !mask cells outside silam domain:
     imask=0
     where ( (center_x < g2%xmax .and. center_x > g2%xmin) .and. (center_y < g2%ymax .and. center_y > g2%ymin) )
          imask=1
     endwhere
     print*,"creando netcdf"
     !create and write netCDF:                                                                              
     call check(nf90_create(gridFile, NF90_CLOBBER, ncid))                                                  
         !! Defino dimensiones:                                                                             
         call check(nf90_def_dim(ncid, "grid_size"     , grid_size    , size_dim_id    ))                   
         call check(nf90_def_dim(ncid, "grid_corners"  , grid_corners , corn_dim_id    ))                   
         call check(nf90_def_dim(ncid, "grid_rank"     , grid_rank    , rank_dim_id    ))                   
         !! Defino variables:                                                                               
         !grid_dims                                                                                         
         call check(nf90_def_var(ncid,"grid_dims"       ,NF90_INT   , [rank_dim_id], var_id             )); 
         !center_lat                                                                                        
         call check(nf90_def_var(ncid,"grid_center_lat" ,NF90_DOUBLE, [size_dim_id], var_id             )); 
         call check(nf90_put_att(ncid, var_id,"units"               , "degrees"       ));                   
         !center_lon                                                                                        
         call check(nf90_def_var(ncid,"grid_center_lon" ,NF90_DOUBLE, [size_dim_id], var_id             )); 
         call check(nf90_put_att(ncid, var_id,"units"               , "degrees"       ));                   
         !imask                                                                                             
         call check(nf90_def_var(ncid,"grid_imask"      ,NF90_INT   , [size_dim_id], var_id             )); 
         call check(nf90_put_att(ncid, var_id,"units"               , "unitless"      ));                   
         !corners_lat                                                                                       
         call check(nf90_def_var(ncid,"grid_corner_lat" ,NF90_DOUBLE, [corn_dim_id,size_dim_id], var_id )); 
         call check(nf90_put_att(ncid, var_id,"units"               , "degrees"       ));                   
         !corners_lon                                                                                       
         call check(nf90_def_var(ncid,"grid_corner_lon" ,NF90_DOUBLE, [corn_dim_id,size_dim_id], var_id  ));
         call check(nf90_put_att(ncid, var_id,"units"               , "degrees"       ));
         !! Defino attributos globales
         call check(nf90_put_att(ncid, nf90_global      ,"title"    , trim(g1%gridName)))
     call check(nf90_enddef(ncid))
     !End NetCDF define mode
     !Abro NetCDF outFile
     call check(nf90_open(gridFile, nf90_write, ncid))
        !grid_dims                                       
        call check(nf90_inq_varid(ncid, 'grid_dims'      , var_id        ))
        call check(  nf90_put_var(ncid, var_id           , [g1%nx,g1%ny] ))
        call check(nf90_inq_varid(ncid, 'grid_center_lat', var_id        ))
        call check(  nf90_put_var(ncid, var_id           , center_y      ))
        call check(nf90_inq_varid(ncid, 'grid_center_lon', var_id        ))
        call check(  nf90_put_var(ncid, var_id           , center_x      ))
        call check(nf90_inq_varid(ncid, 'grid_imask'     , var_id        ))
        call check(  nf90_put_var(ncid, var_id           , imask         ))
        call check(nf90_inq_varid(ncid, 'grid_corner_lat', var_id        ))
        call check(  nf90_put_var(ncid, var_id           , corner_y      ))
        call check(nf90_inq_varid(ncid, 'grid_corner_lon', var_id        ))
        call check(  nf90_put_var(ncid, var_id           , corner_x      ))
     call check(nf90_close(ncid))
     !Cierro NetCDF outFile
  endif

endsubroutine


!***************************!
!* createInterpFile         !
!***************************!
!subroutine createInterpFile(gridFile1,gridFile2,method)
subroutine createInterpFile(g1,g2,method)
   !objetivo: g1,g2, method ==> interp_file = 'remaps/rmp_<g1Name>_to_<g2Name>_<method>.nc'

   !read iFile -> get x,y coords

   !subroutine remap(g,iFile,oFile,varlist, mapMethod)
   implicit none 
   !-----------------------------------------------------------------------
   ! external variables
   !-----------------------------------------------------------------------
   type(regular_grid_type) :: g1,g2   !src & dest grid 
   character(*)    :: method          !bilinear, bicubic, conserv (1st and 2nd).
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
   if (errorCode /= 0) then; print*, 'error initializing SCRIP';stop;endif
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
    mapMethod=method                   ! mapMethod = ’conservative’        !method

    gridFile1="grd_"//trim(g1%gridName)//".nc"
    gridFile2="grd_"//trim(g2%gridName)//".nc"
    interpFile1="rmp_"//trim(g1%gridName)//"_to_"//trim(g2%gridName)//"_"//trim(mapMethod)//".nc"    ! interpFile1 = ’map_1_output_file_name’ !
    interpFile2="rmp_"//trim(g2%gridName)//"_to_"//trim(g1%gridName)//"_"//trim(mapMethod)//".nc"    ! interpFile2 = ’map_2_output_file_name’ ! 
    mapName1=trim(g1%gridName)//"_to_"//trim(g2%gridName)//"_"//trim(mapMethod)                    ! mapName1 = ’name_for_mapping_1’        !
    mapName2=trim(g2%gridName)//"_to_"//trim(g1%gridName)//"_"//trim(mapMethod)                    ! mapName2 = ’name_for_mapping_2’        !

    outputFormat='scrip'               ! outputOpt = ’scrip’               !format of output file. 
                                       
    restrict_type='latitude'           ! restrict_type = ’latitude’        !'latitude' or 'latitude-longitude' search restriction to avoid N^2 search
    num_srch_bins=90                   ! num_srch_bins = 90                !number of bins for search
    normalizeOpt='fracArea'            ! normalizeOpt = ’fracArea’         !'fracArea', 'none', 'dstArea' (normalization for conservative remaping)
    luse_grid1_area=.false.            !
    luse_grid2_area=.false.            !
    npseg=11                           !
    north_thresh=1.5_SCRIP_r8          !
    south_thresh=-2.0_SCRIP_r8         !
    nthreads=2                         !
                                        
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
       print*,'Error: unknown mapping method';stop
    end select

    select case(trim(normalizeOpt))
    case ('none')
       norm_opt = norm_opt_none
    case ('fracArea')
       norm_opt = norm_opt_frcarea
    case ('destArea')
       norm_opt = norm_opt_dstarea
    case default
       print*,'Error: unknown normalization option';stop
    end select
    !-----------------------------------------------------------------------
    ! *create gridFiles (if they don't exist)
    !-----------------------------------------------------------------------
       call grid2gridFile(g1,g2)  !src gridFile 
       call grid2gridFile(g2,g2)  !dst gridFile
    !-----------------------------------------------------------------------
    !  initialize grid information for both grids
    !-----------------------------------------------------------------------
       call grid_init(gridFile1, gridFile2, errorCode)     !src/grids.f
       if (SCRIP_ErrorCheck(errorCode, rtnName, 'Error initializing grids')) then
           call SCRIP_driverExit(errorCode)
       endif
       write(SCRIP_stdout, *) 'Computing remappings between: ',grid1_name
       write(SCRIP_stdout, *) '                         and  ',grid2_name
    !-----------------------------------------------------------------------
    !  initialize some remapping variables.
    !-----------------------------------------------------------------------
       call init_remap_vars                                !src/remap_vars.f
    !-----------------------------------------------------------------------
    !  call appropriate interpolation setup routine based on type of remapping requested.
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
          print*,'Error: Invalid Map Type';stop;
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
    
       call write_remap(mapName1, mapName2, interpFile1, interpFile2, outputFormat, errorCode )
          if (SCRIP_ErrorCheck(errorCode, rtnName, 'error in write_remap')) &
       call SCRIP_driverExit(errorCode)

    !-----------------------------------------------------------------------
    !  All done, exit gracefully
    !-----------------------------------------------------------------------
       call SCRIP_driverExit(errorCode)

end subroutine




!!subroutine applyRemap(var,interp_file)
!subroutine applyRemap(var1,g1,g2,var2,method)
!
!
!      gridFile1 = "grd_"//trim(g1%gridName)//".nc"
!      gridFile2 = "grd_"//trim(g2%gridName)//".nc"
!      interpFile="rmp_"//trim(g1%gridName)//"_to_"//trim(g1%gridName)//"_"//trim(method)//".nc"
!
!      call grid2gridFile(g1,g2)
!      call grid2gridFile(g2,g2)
!      call createInterpFile(g1,g2,method)
!      !-----------------------------------------------------------------------
!      !     initialize SCRIP
!      !-----------------------------------------------------------------------
!      errorCode = SCRIP_Success
!      call SCRIP_CommInitMessageEnvironment
!      call SCRIP_Initialize(errorCode)
!      if (SCRIP_errorCheck(errorCode, rtnName, 'error initializing SCRIP')) then
!        call SCRIP_TestExit(errorCode)
!      endif
!      !-----------------------------------------------------------------------
!      !     read namelist for file and mapping info
!      !-----------------------------------------------------------------------
!      call SCRIP_IOUnitsGet(iunit)
!      open(iunit, file='scrip_test_in', status='old', form='formatted')
!      read(iunit, nml=remap_inputs)
!      call SCRIP_IOUnitsRelease(iunit)
!      write(*,nml=remap_inputs)
!      !-----------------------------------------------------------------------
!      !     read remapping data
!      !-----------------------------------------------------------------------
!       call read_remap(map_name, interp_file, errorCode)
!       if (SCRIP_ErrorCheck(errorCode, rtnName,
!      &    'error reading remap file')) call SCRIP_TestExit(errorCode)
!      !-----------------------------------------------------------------------
!      !     Allocate variables 
!      !-----------------------------------------------------------------------
!      
!      !-----------------------------------------------------------------------
!      !     apply remap        
!      !-----------------------------------------------------------------------
!      if (map_type /= map_type_bicubic) then
!         call remap(grid2_tmp, wts_map1, grid2_add_map1, grid1_add_map1, grid1_array)
!       else
!         call remap(grid2_tmp, wts_map1, grid2_add_map1, grid1_add_map1, grid1_array,
!      &   src_grad1=grad1_lat, src_grad2=grad1_lon, src_grad3=grad1_latlon)
!       endif
!end subroutine





!***********************************************************************
   subroutine SCRIP_driverExit(errorCode)
! !DESCRIPTION:
!  This routine exits the SCRIP driver program. It first calls the 
!  SCRIP error print function to print any errors encountered and then
!  exits the message environment before stopping.
   use SCRIP_KindsMod
   use SCRIP_CommMod
   use SCRIP_ErrorMod
   integer (SCRIP_i4), intent(in) :: &
      errorCode        ! error flag to detect any errors encountered
!-----------------------------------------------------------------------
!  call SCRIP error print function to output any logged errors that
!  were encountered during execution.  Then stop.
!-----------------------------------------------------------------------
   call SCRIP_errorPrint(errorCode, SCRIP_masterTask)
   call SCRIP_CommExitMessageEnvironment
   stop
   end subroutine SCRIP_driverExit
!***********************************************************************

!***********************************************************************
!Just for testing purposes:

!subroutine arrayToNetCDF(g,var)
!   implicit none
!   type(regular_grid_type) :: g
!   real :: var(:,:)
!   real :: lat(:), lon(:)
!   integer :: ncid, row_dim_id,col_dim_id,var_id
!
!   lat=[ (g%ymin+0.5*g%dy+i*g%dy,i=0,ny-1) ]
!   lon=[ (g%xmin+0.5*g%dx+i*g%dx,i=0,nx-1) ] 
!   print*,"creando netcdf"
!   gridFile="test_"//g%gridName//"_var.nc"
!   !create and write netCDF:                                                                              
!   call check(nf90_create(gridFile, NF90_CLOBBER, ncid))                                                  
!       !! Defino dimensiones:           
!       call check(nf90_def_dim(ncid, "lon"      , g%nx, col_dim_id ))
!       call check(nf90_def_dim(ncid, "lat"      , g%ny, row_dim_id ))
!       !! Defino variables:     
!       !lat
!       call check(nf90_def_var(ncid,"lat"         ,NF90_FLOAT, [row_dim_id], var_id          ));
!       call check(nf90_put_att(ncid, pollut_var_id,"units"   , "degrees_north" ));
!       !lon
!       call check(nf90_def_var(ncid,"lon"         ,NF90_FLOAT, [col_dim_id],var_id            ));
!       call check(nf90_put_att(ncid, pollut_var_id,"units"   , "degrees_east"  ));
!       !var   
!       call check(nf90_def_var(ncid, var_list(k)  ,NF90_FLOAT, [col_dim_id,row_dim_id], var_id)) !
!       call check(nf90_put_att(ncid, pollut_var_id,"units"   ,"kg"             ))
!       !! Defino attributos globales
!       call check(nf90_put_att(ncid, nf90_global  ,"title"   , trim(g%gridName)//": var"))
!   call check(nf90_enddef(ncid))
!   !End NetCDF define mode
!   !Abro NetCDF outFile
!   call check(nf90_open(gridFile, nf90_write, ncid))
!      !grid_dims                                       
!      call check(nf90_inq_varid(ncid, 'lon'            , var_id   ))
!      call check(  nf90_put_var(ncid, var_id           , lon      ))
!      call check(nf90_inq_varid(ncid, 'lat'            , var_id   ))
!      call check(  nf90_put_var(ncid, var_id           , lat      ))
!      call check(nf90_inq_varid(ncid, 'var'            , var_id   ))
!      call check(  nf90_put_var(ncid, var_id           , var      ))
!   call check(nf90_close(ncid))
!   !Cierro NetCDF outFile
!end subroutine
!***********************************************************************

end module SCRIP
