program finn2cmaq

  use netcdf
  
  implicit none

  INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
  INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646') 

  type proj_type
     character(16)    :: pName     ! nombre de la proyeccion
     character(7)     :: typ_str   ! Integer code for projection TYPE
     integer          :: typ       ! Integer code for projection TYPE
     integer          :: zone      ! For UTM grids.
     real             :: ref_lat,ref_lon,truelat1,truelat2,stand_lon,pole_lat,pole_lon
     character(125)   :: proj4     ! PROJ4 srs definition.
  end type proj_type

  type grid_type
      character(12)    :: gName        !grid-name
      integer          :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
      real             :: dx,dy        !x-y cell dimension (x_cell, y_cell)
      real             :: xorig,yorig  !lower-left x-y coordinates (in the projected system and units)
      real             :: xmin,ymin,xmax,ymax,xc,yc
      real             :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: i,j,k
  integer :: status,ncid,ncid_out,varid
  logical :: file_exists
  
  character(len=256) :: outfile,inpfile,griddescFile,finnFile,finn_files_dir
  character(len=256) :: command
  character(len=10)  :: chemistry, pollut
  character(len=17), allocatable :: var_list(:)              !lista de polluts
  integer :: nvars

  real, dimension(:,:,:,:), allocatable :: output
  real, dimension(:,:), allocatable :: input
  real, dimension(24) :: diurnal_cycle

  character(len=17) :: start_date_str, end_date_str, current_date_str
  integer ::start_date_s, end_date_s, current_date_s
  character(4) :: YYYY
  character(3) :: DDD
  character(2) :: MM, DD, HH
  integer :: todays_date(8)
  
  !Global parameters:
  !Namelist:
  namelist /parameters/ chemistry,start_date_str,end_date_str,griddescFile,diurnal_cycle,finn_files_dir
  open(7, file='parameters'); read(7,parameters); close(7) !leo namelist

  !GRIDDESC parameters:
  call read_GRIDDESC(griddescFile, proj, grid)

  !Inicializo var_list (polluts) de acuerdo a la quimica seleccionada en namelist
  if ( trim(chemistry) == "GEOSchem"   ) then
    var_list=(/'CO2             ','CO              ','NO              ','NO2             ','SO2             ','NH3             ','CH4             ',&
               'ACET            ','ALD2            ','ALK4            ','BENZ            ','C2H2            ','C2H4            ','C2H6            ',&
               'C3H8            ','CH2O            ','GLYC            ','GLYX            ','HAC             ','MEK             ','MGLY            ',&
               'PRPE            ','TOLU            ','XYLE            ','OC              ','BC              ','PM25            ' /)
  else if ( trim(chemistry) == "MOZ4"   ) then
    var_list=(/'CO2             ','CO              ','NO              ','NO2             ','SO2             ','NH3             ','CH4             ',&
               'VOC             ','ACET            ','ALK1            ','ALK2            ','ALK3            ','ALK4            ','ALK5            ',&
               'ARO1            ','ARO2            ','BALD            ','CCHO            ','CCO_OH          ','ETHENE          ','HCHO            ',&
               'HCN             ','HCOOH           ','HONO            ','ISOPRENE        ','MEK             ','MEOH            ','METHACRO        ',&
               'MGLY            ','MVK             ','OLE1            ','OLE2            ','PHEN            ','PROD2           ','RCHO            ',&
               'TRP1            ','OC              ','BC              ','PM10            ','PM25            '/)
  else if ( trim(chemistry) == "SAPRC99") then
    var_list=(/'CO2             ','CO              ','H2              ','NO              ','NO2             ','SO2             ','NH3             ','CH4             ',&
               'NMOC            ','BIGALD          ','BIGALK          ','BIGENE          ','C10H16          ','C2H4            ','C2H5OH          ','C2H6            ',&
               'C3H6            ','C3H8            ','CH2O            ','CH3CHO          ','CH3CN           ','CH3COCH3        ','CH3COCHO        ','CH3COOH         ',& 
               'CH3OH           ','CRESOL          ','GLYALD          ','HCN             ','HYAC            ','ISOP            ','MACR            ','MEK             ',&
               'MVK             ','TOLUENE         ','HCOOH           ','C2H2            ','OC              ','BC              ','PM10            ','PM25            '/)
   else
        print*, "No existe la química seleccionada: chemistry="//trim(chemistry)//". (Opciones validas: GEOSchem, MOZ4, SAPRC99)"; stop;
   end if
   nvars=size(var_list)

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
    !status=system('mkdir finn_data/')
    inquire(file="finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt", exist=file_exists)
                                                                                                                                                                                      
    if ( .not. file_exists) then
              if ( atoi(YYYY) < todays_date(1)-2 ) then
                        command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//YYYY//"/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
                        status=system(command)
               else
                       command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
                        status=system(command)
               end if
               status=system ("gzip -d finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz")
    end if
    finnFile="finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt"
     
    
    ! Create the NetCDF file                                                        
    call check(nf90_create(outFile, NF90_CLOBBER, ncid))
        ! Define the dimensions
        call check(nf90_def_dim(ncid, "TSTEP"    , 24    , tstep_dimid    )) 
        call check(nf90_def_dim(ncid, "DATE_TIME",  2    , date_time_dimid)) 
        call check(nf90_def_dim(ncid, "COL"      ,grid%nx, col_dimid      )) 
        call check(nf90_def_dim(ncid, "ROW"      ,grid%ny, row_dimid      )) 
        call check(nf90_def_dim(ncid, "LAY"      ,  1    , lay_dimid      )) 
        call check(nf90_def_dim(ncid, "VAR"      ,nvar   , var_dimid      )) 
        ! Define the variables
        do k = 1, nvars             
          pollut=var_list(k) 
          call check(nf90_def_var(ncid, pollut      ,NF90_FLOAT  , [col_dimid, row_dimid], pollut_varid))
          call check(nf90_put_att(ncid, pollut_varid, "long_name", "x-coordinate"))
          call check(nf90_put_att(ncid, pollut_varid, "units"    , "x-coordinate"))
        end do
                                                                                                  
        ! Define the attributes
        call check(nf90_put_att(ncid, "NF_GLOBAL", "IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "EXEC_ID"      , "???????????????? "      ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "FTYPE"        , 1                        ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "SDATE"        , YYYY//DDD                ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "STIME"        , 000000                   ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "WDATE"        , 2023001                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "WTIME"        , 000000                   ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "CDATE"        , 2023001                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "CTIME"        , 000000                   ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "TSTEP"        , 10000                    ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "NTHIK"        , 1                        ))   
        call check(nf90_put_att(ncid, "NF_GLOBAL", "NCOLS"        , grid%nx                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "NROWS"        , grid%ny                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "NLAYS"        , grid%nz                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "NVARS"        , nvars                    ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "GDTYP"        , grid%typ                 ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "P_ALP"        , "-50."                   ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "P_BET"        , "-20."                   ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "P_GAM"        , "-65."                   ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "XCENT"        , proj%ref_lon             ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "YCENT"        , proj%ref_lat             ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "XORIG"        , grid%xorig               ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "YORIG"        , grid%yorig               ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "XCELL"        , grid%dx                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "YCELL"        , grid%dy                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "VGTYP"        , "-9999"                  ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "VGTOP"        , "5000.f"                 ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "VGLVLS"       , "1.f, 0.9938147f"        ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "GDNAM"        , grid%gName               ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "UPNAM"        , "OUTCM3IO"               ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "VAR-LIST"     , var_list                 ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "FILEDESC"     , "Merged emissions output file from Mrggrid" ))
        call check(nf90_put_att(ncid, "NF_GLOBAL", "HISTORY"      , ""                       ))


    ! End the NetCDF define mode
    call check(nf90_enddef(ncid))

    !=================
    !Transformo coordenadas + me quedo con los puntos dentro del dominio.
    !command="cat $finnFile | sed 's/D\([-+]\)/e\1/g;s/,/;/g' | awk -F';' 'NR==1{print $0";wkt"} NR>=2{print $0"POINT("$5*1.0,$4*1.0")"}' > tmp.csv  "
    !call system(command)
    !command="ogr2ogr -f CSV tmp_ll_clipped.csv tmp.csv -clipsrc $lonmin $latmin $lonmax $latmax -lco GEOMETRY='AS_WKT' -lco SEPARATOR='SEMICOLON' "
    !call system(command)
    !command="ogr2ogr -f CSV tmp_xy.csv tmp_ll_clipped.csv -s_srs '$srsInp' -t_srs '$srsOut' -lco GEOMETRY='AS_WKT' -lco SEPARATOR='SEMICOLON' "
    !call system(command)
    !command="sed -i 's///g;s/  */ /g;' tmp_xy.csv "
    !call system(command)
    !=================


    !!Loop over hours of the day
    !do HH = 0, 23

    !   !Memory allocation:
    !   allocate(output(nt,nz,ncol,nrow))  ! Allocate memory for output array
    !   allocate(input(ncol,nrow))         ! Allocate memory for input arrays

    !   ! Open output file for writing
    !   call check(nf90_open(outFile, nf90_write, ncid_out))

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

    !     call check(nf90_inq_varid(ncid_out, pollut, varid))         ! Get output data variable id
    !     call check(nf90_put_var(ncid_out, varid, output(1,1,:,:)))  ! Write output to file
    !   end do
 
    !   !Close output file
    !   call check(nf90_close(ncid_out))     ! Close output file
    !end do

    !!Memory deallocation
    !deallocate(input)    ! Free memory
    !deallocate(output)   ! Free memory

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

 !Funciones para trabajar con DATES
 function date(date_str, fmt_str) result(output)
   implicit none
   character(*), intent(in) :: date_str, fmt_str
   character(50)            :: command
   character(20)            :: output
   command="date -d "//trim(date_str)//" '+"//trim(fmt_str)//"'  > tmp_date.txt"
   call system( trim(command) )
   open(9, file='tmp_date.txt', status='old',action='read'); read(9, '(A)', iostat=status) output;  close(9)
 end function

 !Funciones generales:
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
    character(len=15) :: rtoa
    write(rtoa, '(F15.3)') r
    rtoa = adjustl(rtoa)
 end function

 subroutine read_GRIDDESC(griddescFile, p, g)                                             
        implicit none                                                                          
        character(254),intent(in) :: griddescFile                                                 
        type(proj_type) ,intent(inout) :: p                                                      
        type(grid_type) ,intent(inout) :: g
        character(50)  :: dummyvar
        character(254) :: command        
        
        open(2,file=griddescFile, status='old', action='read')                                     !GRIDDESC:
           read(2,*, iostat=status) dummyvar;                                                      !' '
           read(2,*, iostat=status) p%pName;                                                       !projName
           read(2,*, iostat=status) p%typ,p%truelat1,p%truelat2,p%stand_lon,p%ref_lon,p%ref_lat    !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
           read(2,*, iostat=status) dummyvar;                                                      !' '
           read(2,*, iostat=status) g%gName;                                                       !gridName
           read(2,*, iostat=status) p%pName,g%xorig,g%yorig,g%dx,g%dy,g%ny,g%nx                    !projName xorig yorig xcell ycell nrows ncols
        close(2)

        !Calcular otros parametros:
        if (p%typ == 1 ) then           !Geographic:
                p%typ_str='ll';   p%proj4="+proj=latlong +a=6370000.0 +b=6370000.0"  
        else if ( p%typ == 2 ) then     !Lambert Conformal Conic:
                p%typ_str='lcc';  p%proj4="+proj=lcc +lat_1="//rtoa(p%truelat1)//" +lat_2="//rtoa(p%truelat2)//" +lon_0="//rtoa(p%stand_lon)//" +lat_0="//rtoa(p%ref_lat)//" +a=6370000.0 +b=6370000.0 +units=m"
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
        !Obtener coordenadas del centro de la grilla:
        command="gdaltransform -s_srs 'epsg:4326' -t_srs '"//p%proj4//"' <<< $(echo "//rtoa(p%ref_lon)//" "//rtoa(p%ref_lat)//")  > tmp.txt";
        call system(command)
        open(9, file='tmp_date.txt', status='old',action='read'); read(9, *, iostat=status) g%xc, g%yc, dummyvar;  close(9)

        g%xmin=g%xc - (g%xorig)*(-1);  g%xmax=g%xc + (g%xorig)*(-1)
        g%ymin=g%yc - (g%yorig)*(-1);  g%ymax=g%yc + (g%yorig)*(-1)

        command="gdaltransform -t_srs 'epsg:4326' -s_srs '"//p%proj4//"' <<< $(echo "//rtoa(g%xmin)//" "//rtoa(g%ymin)//") > tmp.txt";
        call system(command)
        open(9, file='tmp_date.txt', status='old',action='read'); read(9, *, iostat=status) g%lonmin, g%latmin, dummyvar;  close(9)
        
        command="gdaltransform -t_srs 'epsg:4326' -s_srs '"//p%proj4//"' <<< $(echo "//rtoa(g%xmax)//" "//rtoa(g%ymax)//") > tmp.txt";
        call system(command)
        open(9, file='tmp_date.txt', status='old',action='read'); read(9, *, iostat=status) g%lonmax, g%latmax, dummyvar;  close(9)

 end subroutine

end program finn2cmaq

!>Namelist
!&parameters
!start_date_str="2019-01-01",	!"%Y-%m-%d
!  end_date_str="2019-01-01",	!"%Y-%m-%d
!chemistry="GEOSchem",
!griddescFile="/home/usuario/runs/papila2019/cmaq/mcip/GRIDDESC",               !path al GRIDDESC
!diurnal_cycle=0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0043,0.0300,0.0600,0.1000,0.1400,0.1700,0.1400,0.1200,0.0900,0.0600,0.0300,0.0043,0.0043,0.0043,0.0043,0.0043
!/

!>Makefile
!.SUFFIXES: .o .f90
!
!FC   = gfortran
!LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
!INC  = -I/usr/include
!FFLAGS = -O2 -ffree-line-length-none
!OBJS = finn2cmaq.o
!
!.f90.o:
!	${FC} ${FFLAGS} -c ${INC} $<
!finn2cmaq: ${OBJS}
!		 ${FC} -o $@ ${OBJS} ${LIBS} 
!clean:
!		rm -f *.o *.mod

