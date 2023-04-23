program prepmegan4cmaq

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

  integer :: status,iostat
  integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,pollut_var_id
  logical :: file_exists


  namelist /control/start_date,end_date,finn_data_directory,griddesc_file,diurnal_cycle







  call date_and_time(values=todays_date)       !fecha de hoy.

  !Leo namelist:
  !open(7, file='parameters'); read(7,parameters); close(7) !leo namelist
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'finn2cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file, proj, grid)    !(!) TO-DO: mejorar esta funcion basado en lo que haga IOAPI

  !Loop over each day
   current_date_s = atoi( date(start_date, "%s") )
       end_date_s = atoi( date(  end_date, "%s") )


!       + CTS `MEGAN_CTS` (*Canopy Type Fractions*) file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains canopy fraction information for six canopy types in one variable, CTS, which is nondimensional and ranges from 0-100. The vegetation types are: needleleaf trees, tropical forest trees, temperate broadleaf trees, shrubs, herbaceous plants, and crops.

!       
!       + LDF `MEGAN_LDF` (*Light Dependence Fractions*) file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains nondimensional light dependence fractions for 4 of the 19 MEGAN chemical species.
!       

!       + EF `MEGAN_EFS`  (emission factors). The MEGAN_EFS file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains emission factors for the 19 MEGAN chemical species.
!       

!       + LAI `MEGAN_LAI` (Leaf Area Index). The MEGAN_LAI file is an I/O API GRDDED3 file that is created using the MEGAN preprocessor. It contains leaf area index that is separate from LAI values used in the rest of CMAQ. By default MEGAN will use this file for LAI, but users can choose to use the LAI values that are read in from MCIP files by setting the environmental variable USE_MEGAN_LAI to N in their run script.


  !Loop por cada día:
  do while (current_date_s <= end_date_s)

    YYYY=date("@"//itoa(current_date_s), "%Y")   !año
      MM=date("@"//itoa(current_date_s), "%m")   !mes
     DDD=date("@"//itoa(current_date_s), "%j")   !dia juliano

    write(*,'(A,A4,A,A3,A,A2,A)') "Day: ",YYYY,"-",DDD," (month: ", MM,")"







    current_date_s=current_date_s + 86400  !siguiente día!
  end do




print*, "========================================="
print*, " prepmegan4cmaq: Completed successfully"
print*, "========================================="

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

end program
