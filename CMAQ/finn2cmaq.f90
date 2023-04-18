program netcdf_copy
  use netcdf
  implicit none

  INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
  INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646') 

  integer :: i,j,k
  integer :: status,ncid,ncid_out,varid
  logical :: file_exists
  
  character(len=256) :: outfile,inpfile,griddescFile,finnFile,finn_files_dir
  character(len=200) :: command
  character(len=10) :: chemistry, pollut
  character(len=17), allocatable :: var_list(:) !lista de polluts
 
  real, dimension(:,:,:,:), allocatable :: output
  real, dimension(:,:), allocatable :: input
  real, dimension(24) :: diurnal_cycle

  integer:: nx,ny,nz,nt,nvars,ncol,nrow
  real :: xorig,yorig,dx,dy

  character(len=17) :: start_date_str, end_date_str, current_date_str
  integer ::start_date_s, end_date_s, current_date_s
  character(4) :: YYYY
  character(3) :: DDD
  character(2) :: MM, DD, HH
  !character(len=10) :: sDate,eDate,sTime,eTime
  integer ::todays_date(8) !start_date(8), end_date(8), current_date(8), 
  !character(4) :: year
  !character(3) :: jday
  !character(2) :: month, day,hour
  

  !Global parameters:
  !Namelist:
  namelist /input_parameters/ chemistry,start_date_str,end_date_str,griddescFile,diurnal_cycle,finn_files_dir
  open(7, file='input_parameters'); read(7,input_parameters); close(7) !leo namelist

  !GRIDDESC parameters: nx, ny, nz, dx, dy, xorig, yorig
  !call get_grid_params_from_griddesc(gridescFile,nx,ny,nz,dx,dy,xorig,yorig)
  nx=197;ny=237;nz=1;dx=20000;dy=20000;     !ncols,nrows,xcel, ycell
  xorig=-1970000;yorig=-2370000;

  !Inicializo var_list (polluts) de acuerdo a la quimica seleccionada en namelist
  if ( trim(chemistry) == "GEOSchem") then
    var_list=(/'CO2             ','CO              ','NO              ','NO2             ','SO2             ','NH3             ','CH4             ',&
               'ACET            ','ALD2            ','ALK4            ','BENZ            ','C2H2            ','C2H4            ','C2H6            ',&
               'C3H8            ','CH2O            ','GLYC            ','GLYX            ','HAC             ','MEK             ','MGLY            ',&
               'PRPE            ','TOLU            ','XYLE            ','OC              ','BC              ','PM25            ' /)
  else if ( trim(chemistry) == "MOZ4" ) then
    var_list=(/'CO2             ','CO              ','NO              ','NO2             ','SO2             ','NH3             ','CH4             ',&
               'VOC             ','ACET            ','ALK1            ','ALK2            ','ALK3            ','ALK4            ','ALK5            ',&
               'ARO1            ','ARO2            ','BALD            ','CCHO            ','CCO_OH          ','ETHENE          ','HCHO            ',&
               'HCN             ','HCOOH           ','HONO            ','ISOPRENE        ','MEK             ','MEOH            ','METHACRO        ',&
               'MGLY            ','MVK             ','OLE1            ','OLE2            ','PHEN            ','PROD2           ','RCHO            ',&
               'TRP1            ','OC              ','BC              ','PM10            ','PM25            '/)
  else if ( trim(chemistry) == "SAPRC99" ) then
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
  !start_date(1) = ichar(start_date_str(1:4)) ; start_date(2) =  ichar(start_date_str(6:7)) ;start_date(3) = ichar(start_date_str(9:10) )
  !  end_date(1) = ichar(  end_date_str(1:4)) ;   end_date(2) =  ichar(  end_date_str(6:7)) ;  end_date(3) = ichar(  end_date_str(9:10) )
  call date_and_time(values=todays_date)       !fecha de hoy.

  !command=char(todays_date(1))//"-"//char(todays_date(2))//"-"//char(todays_date(3))
  !jday=julian_day(command)
  !sdate=julian_day(start_date_str)
  !edate=julian_day(  end_date_str)
  !print*,jday,sdate,edate

  !read(start_date_str, *) start_date
  !read(  end_date_str, *)   end_date
  ! Convert the start date to a Julian day
  !call date_and_time(values=start_date); julian_day = start_date(4) - 1721424.5

  !===========================================================================
  ! Loop over each day
  current_date_s=date_to_seconds(start_date_str)
  end_date_s  = date_to_seconds(end_date_str)   

  do while (current_date_s <= end_date_s)
 
    command=seconds_to_YYYYMMDDD(current_date_s)
    read(command, '(A4, 1X, A2, 1X, A3)') YYYY,MM,DDD !aca falla
    !read('2019 02 123', '(A4, 1X, A2, 1X, A3)') YYYY,MM,DDD  !aca falla

    print*, current_date_s!, YYYY, MM, DD
  
    !!Descargo finn file:
    !!get_finn_data(YYYY,DDD,chemistry,todaysYear,finnFile)
    !!status=system('mkdir finn_data/')
    !!inquire(file="finn_data/GLOB_"//chemistry//"_"//char(YYYY)//char(DDD)//".txt", exist=file_exists)
    !!                                                                                                                                                                                  
    !!if ( .not. file_exists) then
    !!          if ( YYYY < todays_date(1)-2 ) then
    !!                    command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//CHAR(YYYY)//"/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt.gz -P finn_data/"
    !!                    status=system(command)
    !!           else
    !!                   command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt.gz -P finn_data/"
    !!                    status=system(command)
    !!           end if
    !!           status=system ("gzip -d finn_data/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt.gz")
    !!end if
    !!finnFile="finn_data/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt"


    !!Leo finn file:
    !!finnFile="finn_data/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt"
    !open( 1, file=finnFile)
    !!do loop sobre las filas de finnFile
    !close(1) !cierro finnFile
    
    !Filtro puntos dentro del dominio
    !call filter_points_insde_domain()
    
    !Crear nuevo netcdf "outfile", con dimensiones, variables, y atributos.
    
    !! Create the NetCDF file                                                        
    !call check(nf90_create(outFile, NF90_CLOBBER, ncid))
    !! Define the dimensions
    !call check(nf90_def_dim(ncid, "TSTEP"    , 24 , tstep_dimid     ))      
    !call check(nf90_def_dim(ncid, "DATE_TIME",  2 , date_time_dimid ))
    !call check(nf90_def_dim(ncid, "COL"      ,ncol, col_dimid       ))
    !call check(nf90_def_dim(ncid, "ROW"      ,nrow, row_dimid       ))
    !call check(nf90_def_dim(ncid, "LAY"      ,  1 , lay_dimid       ))
    !call check(nf90_def_dim(ncid, "VAR"      ,nvar, var_dimid       ))
    !! Define the variables
    !!do k
    !do k = 1, nvars             
    !  pollut=var_list(k) 
    !  call check(nf90_def_var(ncid, pollut, NF90_FLOAT, [col_dimid, row_dimid], pollut_varid))
    !  call check(nf90_put_att(ncid, pollut_varid, "long_name", "x-coordinate"))
    !  call check(nf90_put_att(ncid, pollut_varid, "units"    , "x-coordinate"))
    !end do
    !                                                                                           
    !! Define the attributes
    !call check(nf90_put_att(ncid, pollut_varid, "long_name", "y-coordinate"))
    !
    !! End the NetCDF define mode
    !call check(nf90_enddef(ncid))


    !!Loop over hours of the day
    !do HH = 0, 23

    !   !Memory allocation:
    !   allocate(output(nt,nz,ncol,nrow))  ! Allocate memory for output array
    !   allocate(input(ncol,nrow))         ! Allocate memory for input arrays

    !   ! Open output file for writing
    !   call check(nf90_open(outfile, nf90_write, ncid_out))

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
     !julian_day = julian_day + 1
  end do
  !===========================================================================
contains

  subroutine check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      write(*,*) nf90_strerror(status)
      stop 'netcdf error'
    end if
  end subroutine check

 !Funciones para trabajar con DATES
 function sys(cmd) result(output)
          implicit none
          character(*), intent(in) :: cmd
          character(200)             :: command
          character(20)              :: output
          command=trim(cmd)//" > tmp.txt"
          call system( trim(command) )
          open(9, file='tmp.txt', status='old',action='read'); read(9, '(A)', iostat=status) output;  close(9)
  end function

  function julian_day(day) result(output)
          implicit none
          character(17) :: day
          character(3) :: output
          output=trim(sys("date -d "//trim(day)//" +%j"))
  end function
 
  function date_to_seconds(day) result(output)
          implicit none
          character(17) :: day
          character(50) :: out_str
          integer :: output
          out_str=trim(sys("date -d "//trim(day)//" '+%s' "))
          read(out_str,'(I50)') output
          !output=ichar(sys("date -d "//trim(day)//" '+%s' "))
  end function

  function seconds_to_YYYYMMDDD(day_sec) result(YYYYMMDDD)
          implicit none
          integer :: day_sec
          character(11) :: YYYYMMDDD
          print*,"date -d @"//char(day_sec)//" '+%Y %m %j'"
          YYYYMMDDD=trim(sys("date -d @"//char(day_sec)//" '+%Y %m %j'"))
          !print*, YYYYMMDDD
  end function                                                             


! subroutine get_finn_data(YYYY,DDD,chemistry,todaysYear,finnFile)
!       implicit none
!       character(len=256), intent(in) :: chemistry
!       integer, intent(in)  :: YYYY, DDD, todaysYear
!       character(len=256), intent(inout) :: finnFile
!       character(len=400) :: command
!       logical :: file_exists
!        
!       status=system('mkdir finn_data/')
!       inquire(file="finn_data/GLOB_"//chemistry//"_"//char(YYYY)//char(DDD)//".txt", exist=file_exists)
!
!       if (file_exists) then
!                
!                   if ( YYYY < todaysYear-2 ) then
!                           command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//CHAR(YYYY)//"/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt.gz -P finn_data/"
!                           status=system(command)
!                  else
!                          command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt.gz -P finn_data/"
!                           status=system(command)
!                   end if
!                                                                                                                                                                     
!                   status=system ("gzip -d finn_data/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt.gz")
!       end if
!       finnFile="finn_data/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt"
! end subroutine

end program netcdf_copy


!>Namelist
!&input_parameters
!start_date_str="2019-01-01"	!"%Y-%m-%d
!  end_date_str="2019-01-01"	!"%Y-%m-%d
!chemistry="GEOSchem"
!griddescFile=""                    !path al GRIDDESC
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


