program finn2cmaq

  use netcdf
  
  implicit none

  INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
  INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646') 

  type proj_type
     character(16)    :: pName     ! nombre de la proyeccion
     character(7)     :: typ_str   ! String  code for projection TYPE
     integer          :: typ       ! Integer code for projection TYPE
     real             :: ref_lat,ref_lon,truelat1,truelat2,stand_lon,pole_lat,pole_lon
     character(125)   :: proj4     ! PROJ4 srs definition.
  end type proj_type

  type grid_type
      character(12)    :: gName        !grid-name
      integer          :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
      real             :: dx,dy        !x-y cell dimension (x_cell, y_cell)
      real             :: xmin,ymin,xmax,ymax,xc,yc
      real             :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status
  integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,pollut_var_id
  logical :: file_exists
  
  character(256) :: command
  character(256) :: outfile,inpfile,griddescFile,finnFile,finn_files_dir
  character(16)  :: chemistry, pollut
  character(17), allocatable :: var_list(:),var_units(:) !lista de polluts
  integer :: nvars
  !finnFile:
  character(len=512) :: header, line
  character(10) :: hvars(6) !day,time,genveg,lati,longi,area
  integer            :: day,time,genveg
  real               :: lati,longi,area
  real, allocatable  :: emis(:)               !valores de emision de cada fila de finnFile.
 
  integer :: i,j,k,ii,ij
  real    :: xi,yi
  real, allocatable  :: data(:,:,:)             !buffer donde meter la grilla con valores de emision [nx,ny,nvars]

  real, dimension(:,:,:,:), allocatable :: output
  real, dimension(:,:), allocatable :: input
  real, dimension(24) :: diurnal_cycle

  character(len=17) :: start_date_str, end_date_str !, current_date_str
  integer :: end_date_s, current_date_s ! start_date_s, 
  character(4) :: YYYY
  character(3) :: DDD
  character(2) :: MM!,DD 
  integer      ::  HH
  integer      :: todays_date(8)
  
  !Global parameters:
  !Namelist:
  namelist /parameters/ chemistry,start_date_str,end_date_str,griddescFile,diurnal_cycle,finn_files_dir
  open(7, file='parameters'); read(7,parameters); close(7) !leo namelist

  !GRIDDESC parameters:
  call read_GRIDDESC(griddescFile, proj, grid)

  !Cuestiones de fecha:
  call date_and_time(values=todays_date)       !fecha de hoy.

  !===========================================================================
  ! Loop over each day
   current_date_s = atoi( date(start_date_str, "%s") )
       end_date_s = atoi( date(  end_date_str, "%s") )

  do while (current_date_s <= end_date_s)
 
    YYYY=date("@"//itoa(current_date_s), "%Y")   !año
      MM=date("@"//itoa(current_date_s), "%m")   !mes
     DDD=date("@"//itoa(current_date_s), "%j")   !dia juliano

    WRITE(*,'(A,A4,A,A3,A,A2,A)') "Dia: ",YYYY,"-",DDD," (mes: ", MM,")"

    !Descargo finn file:
    !call system('mkdir finn_data/')
    inquire(file="finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt", exist=file_exists)
    if ( .not. file_exists) then
              if ( atoi(YYYY) < todays_date(1)-2 ) then
                        command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//YYYY//"/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
                       call system(command)
               else
                       command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
                       call system(command)
               end if
               call system ("gzip -d finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz")
    end if
    finnFile="finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt"
    
    !================================================================
    !Abro finnFile:
    open(1,file=finnFile, status='old', action='read')
                                                                                                                             
    !leo header:                                                          !FinnFile header:
    read(1,'(A)') header                                                  !DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,...,PM25
    nvars = COUNT((/ (header(i:i) == ',', i=1,len(header)) /)) - 6 + 1
    allocate(var_list(nvars))
    allocate(var_units(nvars))
    allocate(emis(nvars))
    read(header,*) hvars,var_list
    
    !Allocato e inicializo grilla con emisiones
    allocate(data(grid%nx,grid%ny,nvars))
    data=0.0

    do  !loop por cada fila de finnFile:
        read(1,*) day,time,genveg,lati,longi,area,emis
     
        if ( (lati > grid%latmax .or. lati < grid%latmin) .or. (longi > grid%lonmax .or. longi < grid%lonmin)  ) then
                continue
        else    
                
                call gdalTransform(lati,longi,xi,yi,'epsg:4326',proj%proj4) !(1) pasar lati y longi a proyectada xi, yi

                !================================================
                ! (!) ESTO FALLA:
                ii=abs(grid%xmax-xi)/abs(grid%xmax-grid%xmin)*grid%nx+grid%nx/2 !posición en la grilla (revisar cuentita) 
                ij=abs(grid%ymax-yi)/abs(grid%ymax-grid%ymin)*grid%ny+grid%ny/2 !posición en la grilla (revisar cuentita) 

                print*,ii,ij,"max dims: ",grid%nx, grid%ny
                !=================================================

                do k=1,nvars
                        pollut=var_list(k)
                        if ( trim(pollut) == "OC" .or. trim(pollut)  == "BC" .or.  trim(pollut) == "PM25" .or. trim(pollut) == "PM10" ) then
                                emis(k) = emis(k) / 3600000.0       !  kg/day -> g/s
                                var_units(k)="g/s"             
                        else
                                emis(k) = emis(k) / 3600.0          !mole/day -> mole/s
                                var_units(k) = "mole/s"           
                        endif
                        data(ii,ij,k) = data(ii,ij,k) + emis(k)
                enddo
        endif
    enddo
    
    !Abro NetCDF
    call check(nf90_open(outFile, nf90_write, ncid))
    do hh=1,24
          ! Create the NetCDF file                                                                                          
          call check(nf90_create(outFile, NF90_CLOBBER, ncid))
              ! Define the dimensions
              call check(nf90_def_dim(ncid, "TSTEP"    ,  1    , tstep_dim_id    )) 
              call check(nf90_def_dim(ncid, "DATE_TIME",  2    , date_time_dim_id)) 
              call check(nf90_def_dim(ncid, "COL"      ,grid%nx, col_dim_id      )) 
              call check(nf90_def_dim(ncid, "ROW"      ,grid%ny, row_dim_id      )) 
              call check(nf90_def_dim(ncid, "LAY"      ,  1    , lay_dim_id      )) 
              call check(nf90_def_dim(ncid, "VAR"      ,nvars  , var_dim_id      )) 
              !! Define the attributes
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "EXEC_ID"      , "???????????????? "      ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "FTYPE"        , 1                        ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "SDATE"        , YYYY//DDD                ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "STIME"        , 000000                   ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "WDATE"        , 2023001                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "WTIME"        , 000000                   ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "CDATE"        , 2023001                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "CTIME"        , 000000                   ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "TSTEP"        , 10000                    ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "NTHIK"        , 1                        ))   
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "NCOLS"        , grid%nx                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "NROWS"        , grid%ny                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "NLAYS"        , grid%nz                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "NVARS"        , nvars                    ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "GDTYP"        , grid%typ                 ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "P_ALP"        , "-50."                   ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "P_BET"        , "-20."                   ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "P_GAM"        , "-65."                   ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "XCENT"        , proj%ref_lon             ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "YCENT"        , proj%ref_lat             ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "XORIG"        , grid%xmin                ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "YORIG"        , grid%ymin                ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "XCELL"        , grid%dx                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "YCELL"        , grid%dy                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "VGTYP"        , "-9999"                  ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "VGTOP"        , "5000.f"                 ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "VGLVLS"       , "1.f, 0.9938147f"        ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "GDNAM"        , grid%gName               ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "UPNAM"        , "OUTCM3IO"               ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "VAR-LIST"     , var_list                 ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "FILEDESC"     , "Merged emissions output file from Mrggrid" ))
              !call check(nf90_put_att(ncid, "NF_GLOBAL", "HISTORY"      , ""                       ))
          ! Define the variables
          do k=1, nvars             
            pollut=var_list(k) 
            call check(nf90_def_var(ncid,trim(pollut)  ,NF90_FLOAT  , [col_dim_id, row_dim_id], pollut_var_id))
            call check(nf90_put_att(ncid, pollut_var_id, "units"    , trim(var_units(k)) ))
          end do
          call check(nf90_enddef(ncid))
          !End the NetCDF define mode   
          !!Open NetCDF outFile
          !call check(nf90_open(outFile, nf90_write, ncid))
          do k=1,nvars
              pollut=var_list(k)
              call check(nf90_inq_varid(ncid, trim(pollut), pollut_var_id))                  ! Get output data variable id
              call check(nf90_put_var(ncid, pollut_var_id, data(:,:,k)*diurnal_cycle(hh) ))  ! Write output to file
          enddo
          
          call check(nf90_close(ncid))
          !Close NetCDF outFile
    end do
    
    deallocate(data)                   !Libero memoria
    deallocate(var_list)               !Libero memoria
    deallocate(var_units)              !Libero memoria
    deallocate(emis)                   !Libero memoria
   
    close(1) !Cierro finnFile
    !================================================================
    
    current_date_s=current_date_s + 86400 
  end do

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

 subroutine read_GRIDDESC(griddescFile, p, g)                                             
        implicit none                                                                          
        character(256),intent(in) :: griddescFile                                                 
        type(proj_type) ,intent(inout) :: p                                                      
        type(grid_type) ,intent(inout) :: g
        character(10) :: dummyvar
        open(2,file=griddescFile, status='old', action='read')                                  !GRIDDESC:
           read(2,*) dummyvar;                                                   !' '
           read(2,*) p%pName;                                                    !projName
           read(2,*) p%typ,p%truelat1,p%truelat2,p%stand_lon,p%ref_lon,p%ref_lat !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
           read(2,*) dummyvar;                                                   !' '
           read(2,*) g%gName;                                                    !gridName
           read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%ny,g%nx                   !projName xorig yorig xcell ycell nrows ncols
        close(2)
        
        !Calcular otros parametros:
        if (p%typ == 1 ) then           !Geographic:
                p%typ_str='ll';   p%proj4="+proj=latlong +a=6370000.0 +b=6370000.0"  
        else if ( p%typ == 2 ) then     !Lambert Conformal Conic:
                p%typ_str='lcc';  p%proj4="+proj=lcc +lat_1="//trim(rtoa(p%truelat1))//" +lat_2="//trim(rtoa(p%truelat2))//" +lon_0="//trim(rtoa(p%stand_lon))//" +lat_0="//trim(rtoa(p%ref_lat))//" +a=6370000.0 +b=6370000.0 +units=m"
        else if ( p%typ == 3 ) then     !General Mercator
                p%typ_str="merc"; p%proj4="+proj=merc +lat_ts=${truelat1} +a=6370000.0 +b=6370000.0"
        else if ( p%typ == 4 ) then     !General tangent Stereografic
                p%typ_str='stere';p%proj4="+proj=stere +lat_0=${ref_lat} +lon_0=${stand_lon} +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
        else if ( p%typ == 5 ) then     !UTM
                p%typ_str='utm';  p%proj4="+proj=utm +zone="  
        else if ( p%typ == 6 ) then     !Polar Secant Stereographic
                p%typ_str='stere';p%proj4="+proj=stere +lat_0=${ref_lat} +lon_0=${stand_lon} +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
        else if ( p%typ == 7 ) then     !Equatorial Mercator
                p%typ_str="merc"; p%proj4="+proj=merc +lat_ts=${truelat1} +a=6370000.0 +b=6370000.0"
        else if ( p%typ == 8 ) then     !Transverse Mercator
                p%typ_str='ll';   p%proj4="+proj=latlong +a=6370000.0 +b=6370000.0"  
        else if ( p%typ == 9 ) then     !Lambert Azimuthal Equal-Area
                print*, "proyección: 9 (Lambert Azimuthal Equal-Area) no soportada por esta aplicación."; stop
        else
                print*, "codigo de proyección invalido.", p%typ; stop
        end if
        
        !Obtener coordenadas del centro de la grilla, min y max:
        g%xc=0.0;g%yc=0.0; g%xmax=(g%xmin)*(-1); g%ymax=(g%ymin)*(-1)

        !transformo boundaries a latlon
        call gdalTransform(g%xmin,g%ymin,g%lonmin,g%latmin,p%proj4,'epsg:4326')
        call gdalTransform(g%xmax,g%ymax,g%lonmax,g%latmax,p%proj4,'epsg:4326')

 end subroutine

 subroutine gdalTransform(x1,y1,x2,y2,srs1,srs2)
        implicit none
        real, intent(in)   :: x1,y1
        real, intent(inout):: x2,y2
        character(*)   :: srs1,srs2
        character(10)  :: ellipsoidh
        character(256) :: command
        command="echo "//rtoa(x1)//" "//rtoa(y1)//" | gdaltransform -s_srs '"//trim(srs1)//"' -t_srs '"//trim(srs2)//"'  > tmp_gdal.txt";
        call system(trim(command))
        !print*,trim(command)
        open(9, file='tmp_gdal.txt', status='old',action='read'); read(9,*, iostat=status) x2, y2, ellipsoidh;  close(9)
 end subroutine

end program finn2cmaq

!>Namelist
!&parameters
!start_date_str="2019-01-01",	!"%Y-%m-%d
!  end_date_str="2019-01-01",	!"%Y-%m-%d
!chemistry="GEOSchem",
!griddescFile="/home/usuario/runs/papila2019/cmaq/mcip/GRIDDESC",               !path al GRIDDESC
!GRIDNAME=
!diurnal_cycle=0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0300,0.0600,0.1000,0.1400,0.1700,0.1400,0.1200,0.0900,0.0600,0.0300,0.0043,0.0043,0.0043,0.0043,0.0043
!/

!>Makefile
!.SUFFIXES: .o .f90
!
!FC   = gfortran
!LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
!INC  = -I/usr/include
!FFLAGS = -O2 -ffree-line-length-none -Wunused
!OBJS = finn2cmaq.o
!
!.f90.o:
!	${FC} ${FFLAGS} -c ${INC} $<
!finn2cmaq: ${OBJS}
!		 ${FC} -o $@ ${OBJS} ${LIBS} 
!clean:
!		rm -f *.o *.mod

!Dependencias: 
!       - NetCDF, 
!       - GDAL/OGR
!       - Shell-tools: date (para descarga: wget, gzip)


!Codigo para reciclar:
!!Loop over hours of the day
!do HH = 0, 23
                                                                                                    
!   !Memory allocation:
!   allocate(output(nt,nz,ncol,nrow))  ! Allocate memory for output array
!   allocate(input(ncol,nrow))         ! Allocate memory for input arrays
                                                                                                    
!   ! Open output file for writing
!   call check(nf90_open(outFile, nf90_write, ncid))
                                                                                                    
!   ! Loop over input files and copy data
!   do k = 1, nvars
!     pollut=trim(var_list(k))
!     inpfile = 'tmp_' // trim(pollut) // '.nc'
!     print*, inpfile
!     !call check(nf90_open( inpfile, nf90_noclobber, ncid))   ! Open input file
!     !call check(nf90_inq_varid(ncid, "Band1", varid))            ! Get input data variable id
!     !call check(nf90_get_var(ncid, varid, input))                ! Read input data
!     output(1,1,:,:) = input                                     ! Copy input data to output array
!     !call check(nf90_close(ncid))                                ! Close input file
                                                                                                    
!     call check(nf90_inq_varid(ncid, pollut, varid))         ! Get output data variable id
!     call check(nf90_put_var(ncid, varid, output(1,1,:,:)))  ! Write output to file
!   end do
                                                                                                    
!   !Close output file
!   call check(nf90_close(ncid))     ! Close output file
!end do
                                                                                                    
!!Memory deallocation
!deallocate(input)    ! Free memory
!deallocate(output)   ! Free memory
