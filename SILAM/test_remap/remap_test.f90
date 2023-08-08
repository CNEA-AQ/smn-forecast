program test

   use netcdf
   use SCRIP
   implicit none
   real, parameter :: pi=3.141593
   real, parameter :: deg2rad=pi/180.0
   type grid_type
     character(12)  :: gName                  !grid-name
     character(100) :: proj4                  !any proj srs string descriptor
     integer        :: nx,ny,nz               !number of cells in x-y direction (ncols, nrows, nlevs)
     real           :: dx,dy                  !x-y cell dimension (x_cell, y_cell)
     real           :: xmin,ymin,xmax,ymax    !bbox and center coorindates (in the proj4 system)
   end type grid_type

   type(grid_type ) :: g
   character(len=256) :: iFile, oFile
   !character(len=256) :: var_list(:)
   character(len=256) :: method

   ifile='test_in.nc'
   ofile='test_ou.nc'

   call makeTestgrid(trim(iFile),1.0,1.0)       !creo archivo de prueba.

   !grids specs:
   g%gName='testgrid'
   g%proj4='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
   g%nx=360
   g%ny=180
   g%dx=1.0
   g%dy=1.0
   g%ymin=-90
   g%xmin=-180
   g%ymax=g%ymin+g%dy*g%ny
   g%xmax=g%xmin+g%dx*g%nx

   !interp. method:
   method='bilinear'


   !call SCRIP_remap(g, iFile, oFile, var_list, method)

   contains


   subroutine makeTestGrid(ncfile,dx,dy)
      implicit none
      character(len=*),intent(in) :: ncfile
      real                :: dx,dy !resolution
      real,allocatable    :: lat(:)
      real,allocatable    :: lon(:)
      real,allocatable    :: var(:,:)
      integer :: i,j,nx,ny
      integer :: ncid, x_dim_id,y_dim_id,var_id

      nx=INT(360/dx)
      ny=INT(180/dy)
      print*,"nx,ny =",nx,ny
      allocate(lat(ny))
      allocate(lon(nx))
      allocate(var(nx,ny))

      lon=(/ (i*dx-180. ,i=1, nx,1)  /) !lon=c(-180:180)
      lat=(/ (i*dy- 90. ,i=1, ny,1)  /) !lat=c(-90:90)
      
      ! Calculate the scalar field values using array operations
      do concurrent (i = 1:nx, j = 1:ny)
          var(i, j) = 2.0 + cos(lon(i)*deg2rad)**2 * cos(2.0*lat(j)*deg2rad)
      end do

      call check(nf90_create(ncFile, NF90_CLOBBER, ncid))
          !! Defino dimensiones
          call check(nf90_def_dim(ncid, "lon" ,   nx   ,   x_dim_id    ))
          call check(nf90_def_dim(ncid, "lat" ,   ny   ,   y_dim_id    ))
          !lat
          call check(nf90_def_var(ncid,"lat"         ,NF90_FLOAT     , [y_dim_id], var_id          ));
          call check(nf90_put_att(ncid, var_id,"units"               , "degrees_north" ));
          call check(nf90_put_att(ncid, var_id,"axis"                , "Y"             ));
          call check(nf90_put_att(ncid, var_id,"long_name"           , "latitude"      ));
          call check(nf90_put_att(ncid, var_id,"standard_name"       , "latitude"      ));
          call check(nf90_put_att(ncid, var_id,"_CoordinateAxisType" , "Lat"           ));
          !lon
          call check(nf90_def_var(ncid,"lon"         ,NF90_FLOAT     , [x_dim_id], var_id          ));
          call check(nf90_put_att(ncid, var_id,"units"               , "degrees_east"  ));
          call check(nf90_put_att(ncid, var_id,"axis"                , "X"             ));
          call check(nf90_put_att(ncid, var_id,"long_name"           , "longitude"     ));
          call check(nf90_put_att(ncid, var_id,"standard_name"       , "longitude"     ));
          call check(nf90_put_att(ncid, var_id,"_CoordinateAxisType" , "Lon"           ));
          !var
          call check(nf90_def_var(ncid,"var"         ,NF90_FLOAT      , [x_dim_id, y_dim_id], var_id          ));
          call check(nf90_put_att(ncid, var_id,"units"               , "g/s"           ));
          call check(nf90_put_att(ncid, var_id,"axis"                , "XY"            ));
          call check(nf90_put_att(ncid, var_id,"long_name"           , "var mass flux" ));
          call check(nf90_put_att(ncid, var_id,"standard_name"       , "var"           ));
      call check(nf90_enddef(ncid))
      !End NetCDF define mode

      !Abro NetCDF outFile
      call check(nf90_open(ncFile, nf90_write, ncid))
         !var
         call check(nf90_inq_varid(ncid, 'var'     , var_id))  !Obtengo id de variable
         call check(nf90_put_var(ncid, var_id, var(:,:))    )  !Escribo valores en NetCDF
         !lat
         call check(nf90_inq_varid(ncid, "lat"     , var_id))
         call check(nf90_put_var(ncid, var_id, lat(:)     ))
         !lon
         call check(nf90_inq_varid(ncid, "lon"     , var_id))
         call check(nf90_put_var(ncid, var_id, lon(:)     ))
      call check(nf90_close(ncid))
      !Cierro NetCDF outFile

   endsubroutine

   subroutine check(status)
     integer, intent(in) :: status
     if (status /= nf90_noerr) then
       write(*,*) nf90_strerror(status)
       stop 'netcdf error'
     end if
   end subroutine check

end program
