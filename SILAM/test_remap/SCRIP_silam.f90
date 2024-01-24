module SCRIP
!
!  SCRIP remaping tool. (Adapted to Silam)
!
  use SCRIP_KindsMod            ! module defining data types
  use SCRIP_CommMod             ! for initializing comm environment
  use SCRIP_ErrorMod            ! SCRIP error checking and logging
  use SCRIP_IOUnitsMod          ! manages I/O units
  use SCRIP_ConfigMod           ! SCRIP configuration module
  use SCRIP_InitMod             ! SCRIP initialization

  use constants                 ! module for common constants
  use timers                    ! CPU timers
  use grids                     ! module with grid information

  use remap_vars                ! common remapping variables
  use remap_conservative        ! routines for conservative remap
  use remap_distance_weight     ! routines for dist-weight remap
  use remap_bilinear            ! routines for bilinear interp
  use remap_bicubic             ! routines for bicubic  interp

  use remap_write               ! routines for remap output
  use remap_read                ! routines for reading remap files
  use remap_mod                 ! module containing remapping routines

  use netcdf 
  use proj_silam   
  
  implicit none
 
  private
  public applyRemap
  public regular_grid_type


  type regular_grid_type
    character(12)  :: gridName             !grid-name
    character(200) :: proj4                !any proj srs string descriptor
    integer        :: nx,ny,nz             !number of cells in x-y direction (ncols, nrows, nlevs)
    real           :: dx,dy                !x-y cell dimension (x_cell, y_cell)
    real           :: xmin,ymin,xmax,ymax  !bbox and center coorindates (in the proj4 system)
  end type regular_grid_type

contains

subroutine grid2gridFile(g1,g2)
  !this subroutine creates a gridFile (necesary to pass to SCRIP to build a interp_file between two grids).
  !it gets two "grid_types" with them grid and proj parameters, and it builds the gridFile.
  implicit none
  type(regular_grid_type) :: g1   !source grid (regular on arbitrary coordinate system)
  type(regular_grid_type) :: g2   !destination grid (silam latlon grid)

  character(100)       :: gridFile
  real(8), allocatable :: center_x(:), center_y(:)
  real(8), allocatable :: corner_x(:,:), corner_y(:,:)
  integer, allocatable :: imask(:)
  integer              :: grid_size, grid_rank,grid_corners
  integer              :: i,j,idx,stat
  integer              :: ncid, size_dim_id,rank_dim_id,corn_dim_id,var_id
  logical              :: file_exist

  gridFile = "grd_"//trim(g1%gridName)//".nc"
  
  print*,'[SCRIP]: Check if grid file'//trim(gridFile)//' exists.'
  inquire(file=trim(gridFile),exist=file_exist)
  if( file_exist ) then; 
     print*,"[SCRIP]: File "//trim(gridFile)//" already exist!!";
     continue;
  else
     print*,'[SCRIP]: creating grid file,'//trim(gridFile)
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

     !compute cells corner coordinates: !! MUST BE ORDERED COUNTERCLOCKWISE !!
     corner_x(1,:)=center_x(:)-0.5*g1%dx     !
     corner_x(2,:)=center_x(:)+0.5*g1%dx     ! 4---------3
     corner_x(3,:)=center_x(:)+0.5*g1%dx     ! |         |
     corner_x(4,:)=center_x(:)-0.5*g1%dx     ! |    o    |
     corner_y(1,:)=center_y(:)-0.5*g1%dy     ! |         |   
     corner_y(2,:)=center_y(:)-0.5*g1%dy     ! 1---------2
     corner_y(3,:)=center_y(:)+0.5*g1%dy     !
     corner_y(4,:)=center_y(:)+0.5*g1%dy     !               

     !transform to latlon
     if (g1%proj4 /= g2%proj4) then
        call proj_trans(g1%proj4, g2%proj4, center_x     , center_y     , grid_size)     !  4---------------3 
        call proj_trans(g1%proj4, g2%proj4, corner_x(1,:), corner_y(1,:), grid_size)     !   \             /
        call proj_trans(g1%proj4, g2%proj4, corner_x(2,:), corner_y(2,:), grid_size)     !    \     o     /
        call proj_trans(g1%proj4, g2%proj4, corner_x(3,:), corner_y(3,:), grid_size)     !     \         /
        call proj_trans(g1%proj4, g2%proj4, corner_x(4,:), corner_y(4,:), grid_size)     !      1-------2
     endif

     !mask cells outside silam domain:
     !imask=1
     imask=0
     where ( (center_x < g2%xmax .and. center_x > g2%xmin) .and. (center_y < g2%ymax .and. center_y > g2%ymin) )
          imask=1
     endwhere
     print*,"[SCRIP] Writing NetCDF with grid data:",trim(gridFile)
     !create and write netCDF:                                                                              
     stat=nf90_create(gridFile, NF90_CLOBBER, ncid)                                                  
         !! Defino dimensiones:                                                                             
         stat=nf90_def_dim(ncid, "grid_size"     , grid_size    , size_dim_id    )                   
         stat=nf90_def_dim(ncid, "grid_corners"  , grid_corners , corn_dim_id    )                   
         stat=nf90_def_dim(ncid, "grid_rank"     , grid_rank    , rank_dim_id    )                   
         !! Defino variables:                                                                               
         !grid_dims                                                                                         
         stat=nf90_def_var(ncid,"grid_dims"       ,NF90_INT   , [rank_dim_id], var_id             ) 
         !center_lat                                                                                      
         stat=nf90_def_var(ncid,"grid_center_lat" ,NF90_DOUBLE, [size_dim_id], var_id             ) 
         stat=nf90_put_att(ncid, var_id,"units"               , "degrees"       )                 
         !center_lon                                                                                      
         stat=nf90_def_var(ncid,"grid_center_lon" ,NF90_DOUBLE, [size_dim_id], var_id             ) 
         stat=nf90_put_att(ncid, var_id,"units"               , "degrees"       )                 
         !imask                                                                                           
         stat=nf90_def_var(ncid,"grid_imask"      ,NF90_INT   , [size_dim_id], var_id             ) 
         stat=nf90_put_att(ncid, var_id,"units"               , "unitless"      )                 
         !corners_lat                                                                                     
         stat=nf90_def_var(ncid,"grid_corner_lat" ,NF90_DOUBLE, [corn_dim_id,size_dim_id], var_id ) 
         stat=nf90_put_att(ncid, var_id,"units"               , "degrees"       )                   
         !corners_lon                                                                                       
         stat=nf90_def_var(ncid,"grid_corner_lon" ,NF90_DOUBLE, [corn_dim_id,size_dim_id], var_id )
         stat=nf90_put_att(ncid, var_id,"units"               , "degrees"       )
         !! Defino attributos globales
         stat=nf90_put_att(ncid, nf90_global      ,"title"    , trim(g1%gridName))
     stat=nf90_enddef(ncid)
     !End NetCDF define mode
     !Abro NetCDF outFile
     stat=nf90_open(gridFile, nf90_write, ncid)
        !grid_dims                                       
        stat=nf90_inq_varid(ncid, 'grid_dims'      , var_id); stat=nf90_put_var(ncid, var_id, [g1%nx,g1%ny] )
        stat=nf90_inq_varid(ncid, 'grid_center_lat', var_id); stat=nf90_put_var(ncid, var_id, center_y      )
        stat=nf90_inq_varid(ncid, 'grid_center_lon', var_id); stat=nf90_put_var(ncid, var_id, center_x      )
        stat=nf90_inq_varid(ncid, 'grid_imask'     , var_id); stat=nf90_put_var(ncid, var_id, imask         )
        stat=nf90_inq_varid(ncid, 'grid_corner_lat', var_id); stat=nf90_put_var(ncid, var_id, corner_y      )
        stat=nf90_inq_varid(ncid, 'grid_corner_lon', var_id); stat=nf90_put_var(ncid, var_id, corner_x      )
     stat=nf90_close(ncid)
     !Cierro NetCDF outFile
     deallocate(center_x,center_y,corner_x,corner_y,imask) !free space OK
  endif

end subroutine grid2gridFile


!=======================================================================================
!subroutine remap_field(var1,var2,g1,g2,method)
subroutine applyRemap(var1,var2,g1,g2,method)
      implicit none
      !input variables:
      real                   , intent(in)    :: var1(:,:) !src grid
      real                   , intent(inout) :: var2(:,:) !dst grid
      type(regular_grid_type), intent(in)    :: g1, g2
      character(*)           , intent(in)    :: method     !bilinear, bicubic, conservative 
      !intermediate variables:
      real (SCRIP_r8), dimension(:), allocatable :: &
                                    grid1_array,    &
                                    grid1_tmp,      &
                                    grad1_lat,      &
                                    grad1_lon,      &
                                    grad1_latlon,   &
                                    grad1_lat_zero, &
                                    grad1_lon_zero, &
                                    grid2_array,    &
                                    grid2_err,      &
                                    grid2_tmp
      integer (SCRIP_i4), dimension(:), allocatable :: grid1_imask, grid2_imask, grid2_count
      !
      character (SCRIP_charLength)     ::  &
               gridFile1   , gridFile2   , &! filename of grid file containing grid2
               interpFile1 , interpFile2 , &! filename for output remap data (map2)
               mapName1    , mapName2    , &! name for mapping from grid2 to grid1
               mapMethod   ,               &! choice for mapping method
               normalizeOpt,               &! option for normalizing weights
               outputFormat                 ! option for output conventions
      character (12), parameter :: rtnName = 'SCRIP_driver'
      integer (SCRIP_i4)        :: errorCode!,iunit !,n ! dummy counter & unit number for input configuration file
      logical                   :: file_exist
      integer :: rss  !memory check


   print*," "
   print*,"[SCRIP]: Asked to map field from grid '",g1%gridName,"' to grid '",g2%gridName,"' using '",trim(method),"' method."
   print*," "
   !-----------------------------------------------------------------------
   !  set-up SCRIP run parameters values ( input namelist )
   !-----------------------------------------------------------------------
   print*,"[SCRIP]: Seting-up remaping parameters.."
   num_maps=1                 !(1: forwards, 2: both directions)
   mapMethod=method           !mapMethod = ’conservative’

   gridFile1="grd_"//trim(g1%gridName)//".nc"
   gridFile2="grd_"//trim(g2%gridName)//".nc"

   mapName1=trim(g1%gridName)//"_to_"//trim(g2%gridName)//"_"//trim(mapMethod)                 
   mapName2=trim(g2%gridName)//"_to_"//trim(g1%gridName)//"_"//trim(mapMethod)                

   interpFile1="rmp_"//trim(mapName1)//".nc" 
   interpFile2="rmp_"//trim(mapName2)//".nc"

   outputFormat='scrip'       !format of output file. 
   restrict_type='latitude'   !'latitude' or 'latlon' search restriction to avoid N^2 search
   num_srch_bins=90           !number of bins for search
   normalizeOpt='fracArea'    !'fracArea', 'none', 'destArea' (normalization for conservative remaping)
   luse_grid1_area=.false.    !
   luse_grid2_area=.false.    !
   npseg=11                   !
   north_thresh=1.5           !1.5_SCRIP_r8  !
   south_thresh=-2.0          !-2.0_SCRIP_r8 !
   nthreads=1                 !
   !-----------------------------------------------------------------------
   select case(mapMethod)
   case ('conservative')
      map_type = map_type_conserv
      !luse_grid_centers = .false.
      luse_grid_centers = .true.
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
      stop '[SCRIP] Error: unknown mapping method'
   end select

   select case(trim(normalizeOpt))
   case ('none')
      norm_opt = norm_opt_none
   case ('fracArea')
      norm_opt = norm_opt_frcarea
   case ('destArea')
      norm_opt = norm_opt_dstarea
   case default
      stop '[SCRIP] Error: unknown normalization option'
   end select

   print*,"[SCRIP]: Check if remapFile exists."
   inquire(file=trim(interpFile1),exist=file_exist)

   if( file_exist ) then; 
      print*,"[SCRIP]: Remaping file "//trim(interpFile1)//" already exist!!";
      !print*,"[SCRIP]: Reading remaping data.."
      call read_remap(mapName1, interpFile1, errorCode)
      continue;

   else
      print*,"[SCRIP]: Creating "//trim(interpFile1)//" remaping file..";
      !-----------------------------------------------------------------------
      ! *create gridFiles (if they don't exist)
      !-----------------------------------------------------------------------
      call grid2gridFile(g1,g2)  !src gridFile 
      call grid2gridFile(g2,g2)  !dst gridFile
      !-----------------------------------------------------------------------
      !  initialize grid information for both grids
      !-----------------------------------------------------------------------
      print*,"[SCRIP]: Initialize grid parameters and compute areas if needed.."
      call grid_init(gridFile1,gridFile2,errorCode)     !src/grids.f
      !-----------------------------------------------------------------------
      !  initialize some remapping variables.
      !-----------------------------------------------------------------------
      print*, "[SCRIP]: Initialize remapping variables .."
      call init_remap_vars                               !src/remap_vars.f
      !-----------------------------------------------------------------------
      !  call appropriate interpolation setup routine based on type of remapping requested.
      !-----------------------------------------------------------------------
      print*, '[SCRIP]: Computing remappings between: ',trim(grid1_name),' and ',trim(grid2_name),'..'
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
            stop '[SCRIP]: Error: Invalid Map Type'
         end select                                     
      !-----------------------------------------------------------------------
      !  reduce size of remapping arrays and write remapping info to file.
      !-----------------------------------------------------------------------
         print*,"[SCRIP]: Writing remaping weights and data on file.."
         if (num_links_map1 /= max_links_map1) then
             call resize_remap_vars(1, num_links_map1-max_links_map1)
         endif
         if ((num_maps > 1) .and. (num_links_map2 /= max_links_map2)) then
              call resize_remap_vars(2, num_links_map2-max_links_map2)
         endif
      !-----------------------------------------------------------------------     
      !  Write remaping file with weights and everithing needed for remaping later.
      !-----------------------------------------------------------------------
         print*,"[SCRIP]: Writing remaping weights and data on file.."
         call write_remap(mapName1, mapName2, interpFile1, interpFile2, outputFormat, errorCode)
   endif
      !-----------------------------------------------------------------------
      !     Allocate variables 
      !-----------------------------------------------------------------------
      allocate (grid1_array    (grid1_size),  &
                grid1_tmp      (grid1_size),  &
                grad1_lat      (grid1_size),  &
                grad1_lon      (grid1_size),  &
                grad1_lat_zero (grid1_size),  &
                grad1_lon_zero (grid1_size),  &
                grid1_imask    (grid1_size),  &
                grid2_array    (grid2_size),  &
                grid2_err      (grid2_size),  &
                grid2_tmp      (grid2_size),  &
                grid2_imask    (grid2_size),  &
                grid2_count    (grid2_size))

      !where (grid1_mask)
      !  grid1_imask = 1
      !elsewhere
      !  grid1_imask = 0
      !endwhere
      !where (grid2_mask)
      !  grid2_imask = 1
      !elsewhere
      !  grid2_imask = 0
      !endwhere

      grid1_array=reshape(var1,[g1%nx*g1%ny])
      !-----------------------------------------------------------------------
      !     apply remap        
      !-----------------------------------------------------------------------
      print*,"[SCRIP]: Applying ",trim(method)," remaping from '",trim(g1%gridName),"' to '",trim(g2%gridName),"'."
      call remap(grid2_tmp, wts_map1, grid2_add_map1, grid1_add_map1, grid1_array)
      !-----------------------------------------------------------------------
      ! save output on 2d array (var2)
      !-----------------------------------------------------------------------
      var2=reshape(grid2_tmp, [g2%nx,g2%ny])
      !-----------------------------------------------------------------------
      ! Free memory
      !-----------------------------------------------------------------------
     print*,"[SCRIP]: Deallocating used variables.."
     if ( allocated(grid1_dims          )) deallocate(grid1_dims          ) !from grid_init (grid.f) 
     if ( allocated(grid1_center_lat    )) deallocate(grid1_center_lat    ) !from grid_init (grid.f) 
     if ( allocated(grid1_center_lon    )) deallocate(grid1_center_lon    ) !from grid_init (grid.f)
     if ( allocated(grid1_area          )) deallocate(grid1_area          ) !from grid_init (grid.f)
     if ( allocated(grid1_frac          )) deallocate(grid1_frac          ) !from grid_init (grid.f)
     if ( allocated(grid1_mask          )) deallocate(grid1_mask          ) !from grid_init (grid.f)
     if ( allocated(grid1_corner_lat    )) deallocate(grid1_corner_lat    ) !from grid_init (grid.f)
     if ( allocated(grid1_corner_lon    )) deallocate(grid1_corner_lon    ) !from grid_init (grid.f)
     if ( allocated(grid2_dims          )) deallocate(grid2_dims          ) !from grid_init (grid.f)
     if ( allocated(grid2_center_lat    )) deallocate(grid2_center_lat    ) !from grid_init (grid.f)
     if ( allocated(grid2_center_lon    )) deallocate(grid2_center_lon    ) !from grid_init (grid.f)
     if ( allocated(grid2_area          )) deallocate(grid2_area          ) !from grid_init (grid.f)
     if ( allocated(grid2_frac          )) deallocate(grid2_frac          ) !from grid_init (grid.f)
     if ( allocated(grid2_mask          )) deallocate(grid2_mask          ) !from grid_init (grid.f)
     if ( allocated(grid2_corner_lat    )) deallocate(grid2_corner_lat    ) !from grid_init (grid.f)
     if ( allocated(grid2_corner_lon    )) deallocate(grid2_corner_lon    ) !from grid_init (grid.f)
     if ( allocated(special_polar_cell1 )) deallocate(special_polar_cell1 ) !from grid_init (grid.f)
     if ( allocated(special_polar_cell2 )) deallocate(special_polar_cell2 ) !from grid_init (grid.f)
     if ( allocated(grid1_area_in       )) deallocate(grid1_area_in       ) !from grid_init (grid.f)
     if ( allocated(grid2_area_in       )) deallocate(grid2_area_in       ) !from grid_init (grid.f)
     if ( allocated(bin_addr1           )) deallocate(bin_addr1           ) !from grid_init (grid.f)
     if ( allocated(bin_addr2           )) deallocate(bin_addr2           ) !from grid_init (grid.f)
     if ( allocated(bin_lats            )) deallocate(bin_lats            ) !from grid_init (grid.f)
     if ( allocated(bin_lons            )) deallocate(bin_lons            ) !from grid_init (grid.f)

     if ( allocated(grid1_bound_box     )) deallocate(grid1_bound_box     ) !from grid_init (grid.f)
     if ( allocated(grid2_bound_box     )) deallocate(grid2_bound_box     ) !from grid_init (grid.f)
     if ( allocated(grid1_centroid_lat  )) deallocate(grid1_centroid_lat  ) !from grid_init (grid.f)
     if ( allocated(grid2_centroid_lat  )) deallocate(grid2_centroid_lat  ) !from grid_init (grid.f)
     if ( allocated(grid1_centroid_lon  )) deallocate(grid1_centroid_lon  ) !from grid_init (grid.f)
     if ( allocated(grid2_centroid_lon  )) deallocate(grid2_centroid_lon  ) !from grid_init (grid.f)

     if ( allocated(wts_map1         )) deallocate(wts_map1               ) !from init_remap_vars (remap_vars.f)
     if ( allocated(wts_map2         )) deallocate(wts_map2               ) !from init_remap_vars (remap_vars.f)
     if ( allocated(grid1_add_map1   )) deallocate(grid1_add_map1         ) !from init_remap_vars (remap_vars.f)
     if ( allocated(grid2_add_map1   )) deallocate(grid2_add_map1         ) !from init_remap_vars (remap_vars.f)
     if ( allocated(grid1_add_map2   )) deallocate(grid1_add_map2         ) !from init_remap_vars (remap_vars.f)
     if ( allocated(grid2_add_map2   )) deallocate(grid2_add_map2         ) !from init_remap_vars (remap_vars.f)

    return
end subroutine

end module
