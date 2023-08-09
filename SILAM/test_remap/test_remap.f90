program test

   use netcdf
   use SCRIP

   implicit none

   type(regular_grid_type ) :: g1,g2
   character(len=256) :: iFile, oFile
   !character(len=256) :: var_list(:)
   character(len=256) :: method

   ifile='in.nc'
   ofile='ou.nc'

   !src grids specs:
   g1%gridName='testgrid'
   g1%proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" !'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
   g1%nx=360
   g1%ny=180
   g1%dx=1.0
   g1%dy=1.0
   g1%ymin=-90
   g1%xmin=-180
   
   call makeTestgrid(trim(iFile),g1) !dummy test file.

   !dst grid specs (silam grid):
   g2%gridName="silamgrid"
   g2%proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
   g2%nx=50
   g2%ny=63
   g2%xmin=-72.0
   g2%ymin=-54.0
   g2%dx=0.4
   g2%dy=0.4

   !!create gridFiles
   !call grid2gridFile(g1,g2,"grd_"//trim(g1%gridName)//".nc")
   !call grid2gridFile(g2,g2,"grd_"//trim(g2%gridName)//".nc")

   !!interp. method:
   method='bilinear'

   call createInterpFile(g1,g2,method)

   !!apply remaping
   !call applyRemap(var,interp_file)

   contains


   subroutine makeTestGrid(ncfile,g)
      implicit none
      character(len=*),intent(in) :: ncfile
      type(regular_grid_type ) :: g
      !real                :: dx,dy !resolution
      real,allocatable    :: lat(:)
      real,allocatable    :: lon(:)
      real,allocatable    :: var(:,:)
      integer :: i,j!,nx,ny
      integer :: ncid, x_dim_id,y_dim_id,var_id

      allocate(lat(g%ny))
      allocate(lon(g%nx))
      allocate(var(g%nx,g%ny))

      lon=[ (g%xmin+i*g%dx, i=1, g%nx,1) ] !lon=c(-180:180)
      lat=[ (g%ymin+i*g%dy, i=1, g%ny,1) ] !lat=c(-90:90)
      
      ! Calculate the scalar field values using array operations
      do concurrent (i = 1:g%nx, j = 1:g%ny)
          var(i, j) = 2.0 + cos(lon(i)*deg2rad)**2 * cos(2.0*lat(j)*deg2rad)
      end do

      call check(nf90_create(ncFile, NF90_CLOBBER, ncid))
          !! Defino dimensiones
          call check(nf90_def_dim(ncid, "lon" ,   g%nx   ,   x_dim_id    ))
          call check(nf90_def_dim(ncid, "lat" ,   g%ny   ,   y_dim_id    ))
          !lat
          call check(nf90_def_var(ncid,"lat"         ,NF90_FLOAT     , [y_dim_id], var_id          ));
          call check(nf90_put_att(ncid, var_id,"units"               , "degrees_north" ));
          call check(nf90_put_att(ncid, var_id,"long_name"           , "latitude"      ));
          !lon
          call check(nf90_def_var(ncid,"lon"         ,NF90_FLOAT     , [x_dim_id], var_id          ));
          call check(nf90_put_att(ncid, var_id,"units"               , "degrees_east"  ));
          call check(nf90_put_att(ncid, var_id,"long_name"           , "longitude"     ));
          !var
          call check(nf90_def_var(ncid,"var"         ,NF90_FLOAT      , [x_dim_id, y_dim_id], var_id          ));
          call check(nf90_put_att(ncid, var_id,"units"               , "g/s"           ));
          call check(nf90_put_att(ncid, var_id,"long_name"           , "var mass flux" ));
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

   !subroutine check(status)
   !  integer, intent(in) :: status
   !  if (status /= nf90_noerr) then
   !    write(*,*) nf90_strerror(status)
   !    stop 'netcdf error'
   !  end if
   !end subroutine check

end program
