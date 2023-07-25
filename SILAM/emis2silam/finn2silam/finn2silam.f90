program finn2silam

  use netcdf

  implicit none

  INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
  INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646') 

  type grid_type
      character(12)   :: gName        !grid-name
      integer         :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
      real            :: dx,dy        !x-y cell dimension (x_cell, y_cell)
      real            :: xmin,ymin,xmax,ymax,xc,yc
      real            :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(grid_type) :: grid

  integer :: status,iostat
  integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,pollut_var_id
  logical :: file_exists
  character(256) :: command     !(por si uso wget)
  character(256) :: griddesc_file,finn_data_directory,finnFile,outFile
  character(16)  :: gridname, chemistry, pollut
  character(16), allocatable :: var_list(:), var_units(:) !lista de polluts
  character(800) :: var_list_string
  character(80) :: att_var_desc
  integer :: nvars
  
  !finnFile:
  character(len=512) :: header                !
  character(10)      :: colnames(6)           !columnames de finnFile
  integer            :: day,time,genveg       !vars de finnFile
  real               :: lati,longi,area       !vars de finnFile
  real, allocatable  :: emis(:)               !emision de cada fila de finnFile.
 
  integer :: i,j,k,h,ii,ij                    !indices.
  !real    :: xi,yi
  real   , allocatable  :: data(:,:,:,:,:)         !buffer donde meter la grilla con valores de emision [nt,nz,nx,ny,nvars]
  real   , allocatable  :: lat(:),lon(:),height(:) !buffer donde meter las coordenadas [nx] y [ny]
  integer :: timevar(24)                           !buffer donde meter los valores de tiempo [24]

  real, dimension(24) :: diurnal_cycle

  character(len=17) :: start_date, end_date
  integer :: end_date_s, current_date_s
  character(19) :: current_date !current date with format %Y-%m-%d %H:%M:%S
  character(4) :: YYYY
  character(3) :: DDD
  character(2) :: MM,HH,DD 
  integer      :: todays_date(8)

  real    :: lon_start,lat_start,lon_end,lat_end,dx,dy
  integer ::nx,ny

  namelist /control/ chemistry,start_date,end_date,finn_data_directory,lon_start,lat_start,lon_end,lat_end,dx,dy,nx,ny,diurnal_cycle
  
  call date_and_time(values=todays_date)       !fecha de hoy.

  nx=-999;ny=-999;lon_end=-999.;lat_end=-999.
  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'finn2silam: failed to read namelist; error = ',iostat
    stop
  end if

  !Set grid parameters:
  grid%lonmin=lon_start      ; grid%latmin=lat_start
  grid%dx=dx                 ; grid%dy= dy      
  if ( nx == -999 ) then
        if ( lon_end == -999.0 ) then
            print*, "Error: nx nor lon_end has been initialized on namelist!"; stop;
        else 
            grid%lonmax=lon_end; grid%nx=(grid%lonmax-grid%lonmin)/grid%dx
        endif
  else if ( lon_end == -999.0 ) then
        grid%nx=nx; grid%lonmax=grid%lonmin+grid%nx*grid%dx
  endif
  if ( ny == -999 ) then
        if ( lat_end == -999.0 ) then
            print*, "Error: nx nor lon_end has been initialized on namelist!"; stop;
        else 
            grid%latmax=lat_end; grid%ny=(grid%latmax-grid%latmin)/grid%dy
        endif
  else if ( lat_end == -999.0 ) then
        grid%ny=nx; grid%latmax=grid%latmin+grid%ny*grid%dy
  endif

  allocate(height(1))
  allocate(lon(grid%nx))
  allocate(lat(grid%ny))
  height=[50.0] !50 meters as emission height (change!)
  lon=[(grid%lonmin+i*grid%dx,i=1,grid%nx)]
  lat=[(grid%latmin+i*grid%dy,i=1,grid%ny)]

  !Loop over each day
   current_date_s = atoi( date(start_date, "%s") )
       end_date_s = atoi( date(  end_date, "%s") )

  do while (current_date_s <= end_date_s)
                    
    current_date=date("@"//itoa(current_date_s), "%Y-%m-%d %H:%M:%S")

    YYYY=date("@"//itoa(current_date_s), "%Y")   !año
      MM=date("@"//itoa(current_date_s), "%m")   !mes
      DD=date("@"//itoa(current_date_s), "%d")   !dia juliano
     DDD=date("@"//itoa(current_date_s), "%j")   !dia juliano

    write(*,'(A,A4,A,A3,A,A2,A)') "Day: ",YYYY,"-",DDD," (month: ", MM,")"
    
    finnFile=trim(finn_data_directory)//"GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt"
    inquire(file=finnFile, exist=file_exists)

#ifdef WGET   
    !Descargo finn file:
    if ( .not. file_exists) then
      if ( atoi(YYYY) <= todays_date(1)-2 ) then
         command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//YYYY//"/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
         call system(command)
      else
         command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
         call system(command)
      end if
      call system ("gzip -d finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz")
      inquire(file=finnFile, exist=file_exists)
    end if
#endif
    if ( file_exists ) then
       
       print*," Reading Finn file: ",trim(finnFile)
       !Abro finnFile:
       open(1,file=finnFile, status='old', action='read')
          read(1,'(A)') header                                                  !Aca asumo que FinnFile header es siempre:
          nvars = COUNT((/ (header(i:i) == ',', i=1,len(header)) /)) - 6 + 1    !DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,...,PM25

          allocate(var_list(nvars))                         !array con nombre de polluts
          allocate(emis(nvars))                             !array con emisiones de los polluts
          allocate(var_units(nvars))                        !array con unidades de emisiones
          allocate(data(grid%nx,grid%ny,1,24,nvars))        !grilla de emisiones
          !allocate(timevar(24))                       !variable timevar
          data=0.0

          read(header,*) colnames,var_list
          write(var_list_string,*) var_list                !este es un global attr importante.

          iostat=0
          do while(iostat == 0)  !loop por cada fila de finnFile:
              read(1,*,iostat=iostat) day,time,genveg,lati,longi,area,emis
           
              if ( lati > grid%latmax .or. lati < grid%latmin .or. longi > grid%lonmax .or. longi < grid%lonmin ) then
                continue
              else    
                ii=floor((longi-grid%lonmin)/(grid%lonmax-grid%lonmin)*grid%nx) !calculo posición-X en la grilla
                ij=floor(( lati-grid%latmin)/(grid%latmax-grid%latmin)*grid%ny) !calculo posición-Y en la grilla
                do k=1,nvars
                        pollut=var_list(k)
                        if ( trim(pollut) == "OC" .or. trim(pollut)  == "BC" .or.  trim(pollut) == "PM25" .or. trim(pollut) == "PM10" ) then
                                emis(k) = emis(k) / 3600000.0      !kg/day -> g/s      (estrictamente es cierto cuando aplico el ciclo diurno)
                                var_units(k)="g/s"             
                        else
                                emis(k) = emis(k) / 3600.0         !mole/day -> mole/s (estrictamente es cierto cuando aplico el ciclo diurno)
                                var_units(k) = "mole/s"           
                        endif

                        do h=1,24
                            timevar(h) = current_date_s + 3600*(h-1)  !seconds from base_date_s
                            data(ii,ij,1,h,k) = data(ii,ij,1,h,k) + emis(k)*diurnal_cycle(h)
                        enddo
                enddo
              endif
          enddo
       close(1) !Cierro finnFile
       
       outFile="./emis_fires_"//YYYY//MM//DD//"_d01.nc"

       print*," Creating NetCDF file: ",trim(outFile)
       ! Create the NetCDF file                                                                                          
       call check(nf90_create(outFile, NF90_CLOBBER, ncid))
           !! Defino dimensiones
           call check(nf90_def_dim(ncid, "time"    ,   24   , tstep_dim_id    )) !24 because daily file, hourly time-step.
           call check(nf90_def_dim(ncid, "lon"     , grid%nx, col_dim_id      )) 
           call check(nf90_def_dim(ncid, "lat"     , grid%ny, row_dim_id      )) 
           call check(nf90_def_dim(ncid, "height"  ,   1    , lay_dim_id      )) 
           !! Defino attributos globales
           call check(nf90_put_att(ncid, nf90_global,"title"       , "SILAM_OUTPUT"   ))
           call check(nf90_put_att(ncid, nf90_global,"Conventions" , "CF-1.3"         ))
           call check(nf90_put_att(ncid, nf90_global,"source"      ,"SILAM v5_7_SVN (r588062)" ))
           call check(nf90_put_att(ncid, nf90_global,"_CoordinateModelRunDate",current_date//"Z"))
           call check(nf90_put_att(ncid, nf90_global,"grid_projection" , "lonlat"     ))
           call check(nf90_put_att(ncid, nf90_global,"pole_lat", -90.0                ))
           call check(nf90_put_att(ncid, nf90_global,"pole_lon",   0.0                ))
           call check(nf90_put_att(ncid, nf90_global,"history" , ""                   ))
           !!Defino variables extra (height,lat,lon,time)
           call check(nf90_def_var(ncid,"height"      ,NF90_FLOAT            , [lay_dim_id], pollut_var_id          ));
           call check(nf90_put_att(ncid, pollut_var_id, "units"              , "m"                                  ));
           call check(nf90_put_att(ncid, pollut_var_id, "long_name"          , "layer midpoint constant height from surface"));
           call check(nf90_put_att(ncid, pollut_var_id, "axis"               , "Z"                                  ));
           call check(nf90_put_att(ncid, pollut_var_id, "standard_name"      ,"layer_midpoint_height_above_ground"))
           !lat
           call check(nf90_def_var(ncid,"lat"         ,NF90_FLOAT            , [row_dim_id], pollut_var_id          ));
           call check(nf90_put_att(ncid, pollut_var_id,"units"               , "degrees_north" ));
           call check(nf90_put_att(ncid, pollut_var_id,"axis"                , "Y"             ));        
           call check(nf90_put_att(ncid, pollut_var_id,"long_name"           , "latitude"      )); 
           call check(nf90_put_att(ncid, pollut_var_id,"standard_name"       , "latitude"      )); 
           call check(nf90_put_att(ncid, pollut_var_id,"_CoordinateAxisType" , "Lat"           ));
           !lon
           call check(nf90_def_var(ncid,"lon"         ,NF90_FLOAT            , [col_dim_id], pollut_var_id          ));
           call check(nf90_put_att(ncid, pollut_var_id,"units"               , "degrees_east"  ));
           call check(nf90_put_att(ncid, pollut_var_id,"axis"                , "X"             ));        
           call check(nf90_put_att(ncid, pollut_var_id,"long_name"           , "longitude"     )); 
           call check(nf90_put_att(ncid, pollut_var_id,"standard_name"       , "longitude"     )); 
           call check(nf90_put_att(ncid, pollut_var_id,"_CoordinateAxisType" , "Lon"           ));
           !time
           call check(nf90_def_var(ncid,"time"        ,NF90_INT              , [tstep_dim_id], pollut_var_id  ));
           call check(nf90_put_att(ncid, pollut_var_id,"units"               , "seconds since "//current_date//" UTC" ));
           call check(nf90_put_att(ncid, pollut_var_id,"long_name"           , "time"                                  ));        
           call check(nf90_put_att(ncid, pollut_var_id,"axis"                , "T"                                     )); 
           call check(nf90_put_att(ncid, pollut_var_id,"calendar"            , "standard"                              )); 
           call check(nf90_put_att(ncid, pollut_var_id,"standard_name"       , "time"                                  ));

           !emission variables:
           do k=1, nvars             
             pollut=var_list(k) 
             att_var_desc=trim(pollut)//"[1]"
             call check(nf90_def_var(ncid, pollut       ,NF90_FLOAT , [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], pollut_var_id)) !
             call check(nf90_put_att(ncid, pollut_var_id,"units"    , trim(var_units(k)) ))
             call check(nf90_put_att(ncid, pollut_var_id,"long_name", pollut             ))
             call check(nf90_put_att(ncid, pollut_var_id,"var_desc" , att_var_desc       ))
           end do

       call check(nf90_enddef(ncid))
       !End NetCDF define mode   
        
       !Abro NetCDF outFile
       call check(nf90_open(outFile, nf90_write, ncid))
       do k=1,nvars
           pollut=var_list(k)
           call check(nf90_inq_varid(ncid, trim(pollut), pollut_var_id))     !Obtengo id de variable
           call check(nf90_put_var(ncid, pollut_var_id, data(:,:,:,:,k)  ))  !Escribo valores en NetCDF
       enddo
       !height
       call check(nf90_inq_varid(ncid, "height"  , pollut_var_id))                                
       call check(nf90_put_var(ncid, pollut_var_id, height(:)   ))
       !lat
       call check(nf90_inq_varid(ncid, "lat"     , pollut_var_id))                                
       call check(nf90_put_var(ncid, pollut_var_id, lat(:)     ))
       !lon
       call check(nf90_inq_varid(ncid, "lon"     , pollut_var_id))                                
       call check(nf90_put_var(ncid, pollut_var_id, lon(:)     ))
       !time
       call check(nf90_inq_varid(ncid, "time"    , pollut_var_id))                                
       call check(nf90_put_var(ncid, pollut_var_id, timevar(:) ))

       call check(nf90_close(ncid))
       !Cierro NetCDF outFile

    deallocate(data)       !Libero memoria
    timevar=0
    deallocate(var_list)   !Libero memoria
    deallocate(var_units)  !Libero memoria
    deallocate(emis)       !Libero memoria
   
    endif
    
    current_date_s=current_date_s + 86400  !siguiente día!
  end do

print*, "==================================="
print*, " finn2silam: Completed successfully "
print*, "==================================="
contains

 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

 !Interfaz a "date"
 function date(date_str, fmt_str) result(output)
   implicit none
   character(*), intent(in) :: date_str, fmt_str
   character(256)           :: command
   character(20)            :: output
   command="date -d "//trim(date_str)//" '+"//trim(fmt_str)//"'  > tmp_date.txt"
   call system( trim(command) )
   !print*,trim(command)
   open(9, file='tmp_date.txt', status='old',action='read'); read(9, '(A)', iostat=status) output;  close(9)
   call system('rm tmp_date.txt')
 end function

 function atoi(str)     !string -> int
   implicit none
   character(len=*), intent(in) :: str
   integer :: atoi
   read(str,*) atoi
 end function
 function itoa(i)       !int -> string
    implicit none
    integer, intent(in) :: i
    character(len=20) :: itoa
    write(itoa, '(i0)') i
    itoa = adjustl(itoa)
 end function
 function rtoa(r)       !real -> string
    implicit none
    real, intent(in) :: r
    character(len=16) :: rtoa
    write(rtoa, '(F16.3)') r
    rtoa = adjustl(rtoa)
 end function

end program finn2silam
