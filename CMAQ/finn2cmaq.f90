program netcdf_copy
  use netcdf
  implicit none

  integer :: status,ncid,ncid_out,varid,nrow,ncol,nt,nz,i,j,k,nvars
  integer, dimension(:,:), allocatable :: data
  character(len=16) :: var_list(50)
  character(len=12) :: pollut, chemistry
  character(len=256) :: outfile, inpfile,gridescFile,finnFile
  real, dimension(:,:,:,:), allocatable :: output
  real, dimension(:,:), allocatable :: input
  real, dimension(23) :: diurnal_cycle
  character(len=19) :: start_date_str, end_date_str, current_date_str
  integer :: start_date(8), end_date(8), current_date(8), todays_date(8)
  integer :: YYYY, MM, DD, HH, DDD
  real :: julian_day

  !Global parameters:
  !Namelist:
  namelist /input_parameters/ chemistry,start_date_str,end_date_str,gridescFile,diurnal_cycle
  open(7, file='input_parameters'); read(7, nml=input_parameters); close(7) !leo namelist

  !GRIDDESC parameters: nx, ny, nz, dx, dy, xorig, yorig
  !call get_grid_params_from_griddesc(gridescFile,nx,ny,nz,dx,dy,xorig,yorig)
  nx=197;ny=237;nz=1;dx=20000;dy=20000;     !ncols,nrows,xcel, ycell
  xorig=-1970000;yorig=-2370000;

  !Inicializo var_list (polluts) de acuerdo a la quimica seleccionada en namelist
  call get_var_list(chemistry,var_list,nvars)

  !Cuestiones de fecha:
  call date_and_time(values=todays_date)       !fecha de hoy.
  read(start_date_str, *) start_date
  read(end_date_str, *) end_date
  ! Convert the start date to a Julian day
  call date_and_time(values=start_date); julian_day = start_date(4) - 1721424.5

  !===========================================================================
  ! Loop over each day
  do while (julian_day < end_date(4) - 1721424.5)
 
    call julian_to_date(julian_day,YYYY,MM,DD)
  
    !!Descargo finn file:
    !!get_finn_data(YYYY,DDD,chemistry,todaysYear,finnFile)

    !!Leo finn file:
    !!finnFile="finn_data/GLOB_"//chemistry//"_"//CHAR(YYYY)//CHAR(DDD)//".txt"
    !open( 1, file=finnFile)
    !!do loop sobre las filas de finnFile
    !close(1) !cierro finnFile
    
    !Filtro puntos dentro del dominio
    !call filter_points_insde_domain()
    
    do HH = 0, 23
   
        !Crear nuevo netcdf "outfile", con dimensiones, variables, y atributos.

        !Memory allocation:
        allocate(output(nt,nz,ncol,nrow))  ! Allocate memory for output array
        allocate(input(ncol,nrow))         ! Allocate memory for input arrays

        ! Open output file for writing
        call check(nf90_open(outfile, nf90_write, ncid_out))

        ! Loop over input files and copy data
        do k = 1, nvars
          pollut=var_list(k)
          inpfile = 'tmp_' // trim(pollut) // '.nc'
          print*, inpfile
          !call check(nf90_open( inpfile, nf90_noclobber, ncid))   ! Open input file
          !call check(nf90_inq_varid(ncid, "Band1", varid))            ! Get input data variable id
          !call check(nf90_get_var(ncid, varid, input))                ! Read input data
          output(1,1,:,:) = input                                     ! Copy input data to output array
          !call check(nf90_close(ncid))                                ! Close input file

          call check(nf90_inq_varid(ncid_out, pollut, varid))         ! Get output data variable id
          call check(nf90_put_var(ncid_out, varid, output(1,1,:,:)))  ! Write output to file
        end do
 
        !Close output file
        call check(nf90_close(ncid_out))     ! Close output file
     end do

     !Memory deallocation
     deallocate(input)    ! Free memory
     deallocate(output)   ! Free memory

     julian_day = julian_day + 1
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

 subroutine get_var_list(chemistry,var_list,nvars)
        implicit none
        character, intent(in) :: chemistry
        character(len=16), intent(inout) :: var_list(50)
        integer, intent(inout) :: nvars
        
        if ( chemistry == "GEOSchem") then
              var_list=(/ 'CO2','CO','NO','NO2','SO2','NH3','CH4','ACET','ALD2','ALK4','BENZ','C2H2','C2H4','C2H6','C3H8','CH2O','GLYC','GLYX','HAC','MEK','MGLY','PRPE','TOLU','XYLE','OC','BC','PM25' /)
              nvars=27
        else if ( chemistry == "MOZ4" ) then
              var_list=(/ 'CO2','CO','NO','NO2','SO2','NH3','CH4','VOC','ACET','ALK1','ALK2','ALK3','ALK4','ALK5','ARO1','ARO2','BALD','CCHO','CCO_OH','ETHENE','HCHO','HCN','HCOOH','HONO','ISOPRENE','MEK','MEOH','METHACRO','MGLY','MVK','OLE1','OLE2','PHEN','PROD2','RCHO','TRP1','OC','BC','PM10','PM25' /)
              nvars=40
        else if ( chemistry == "SAPRC99" ) then
              var_list=(/ 'CO2','CO','H2','NO','NO2','SO2','NH3','CH4','NMOC','BIGALD','BIGALK','BIGENE','C10H16','C2H4','C2H5OH','C2H6','C3H6','C3H8','CH2O','CH3CHO','CH3CN','CH3COCH3','CH3COCHO','CH3COOH','CH3OH','CRESOL','GLYALD','HCN','HYAC','ISOP','MACR','MEK','MVK','TOLUENE','HCOOH','C2H2','OC','BC','PM10','PM25' /)
              nvars=40
         else
              print*, "No existe la química seleccionada: chemistry="//chemistry//". (Opciones validas: GEOSchem, MOZ4, SAPRC99)"; stop;
         end if
 end subroutine get_var_list


 subroutine create_new_netcdf(outFile, ncol, nrow, nz, nt, var_list,start_date)
        implicit none
        character(len=256), intent(in) :: outFile
        integer, intent(in) :: ncol,nrow, nvar
        integer ncid,tstep_dimid,date_time_dimid,col_dimid,row_dimid,lay_dimid,var_dimid, pollut_varid
        ! Create the NetCDF file                                                        
        call check(nf90_create(outFile, NF90_CLOBBER, ncid))
        
        ! Define the dimensions
        call check(nf90_def_dim(ncid, "TSTEP"    , 24 , tstep_dimid     ))      
        call check(nf90_def_dim(ncid, "DATE_TIME",  2 , date_time_dimid ))
        call check(nf90_def_dim(ncid, "COL"      ,ncol, col_dimid       ))
        call check(nf90_def_dim(ncid, "ROW"      ,nrow, row_dimid       ))
        call check(nf90_def_dim(ncid, "LAY"      ,  1 , lay_dimid       ))
        call check(nf90_def_dim(ncid, "VAR"      ,nvar, var_dimid       ))
        
        ! Define the variables
        !do k
        do k = 1, nvars             
          pollut=var_list(k) 
          call check(nf90_def_var(ncid, pollut, NF90_FLOAT, [col_dimid, row_dimid], pollut_varid))
          call check(nf90_put_att(ncid, pollut_varid, "long_name", "x-coordinate"))
          call check(nf90_put_att(ncid, pollut_varid, "units"    , "x-coordinate"))
        end do

        ! Define the attributes
        call check(nf90_put_att(ncid, pollut_varid, "long_name", "y-coordinate"))
        
        ! End the NetCDF define mode
        call check(nf90_enddef(ncid))

 end subroutine



end program netcdf_copy


!>Namelist
!&input_parameters
!start_date="2019-01-01"	!"%Y-%m-%d %H"
!  end_date="2019-01-01"	!"%Y-%m-%d %H"
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
!		 ${F90} -o $@ ${OBJS} ${LIBS} 
!clean:
!		rm -f *.o *.mod







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

