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






  namelist /control/start_date,end_date,griddesc_file
  namelist /inp_files/crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file, nl_tree_frac_file, bl_tree_frac_file, tp_tree_frac_file,ecotype_file,lai_file,veg_cov_file,wrfout_file

  call date_and_time(values=todays_date)       !fecha de hoy.

  !Leo namelist:
  !open(7, file='parameters'); read(7,parameters); close(7) !leo namelist
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
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
           read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny                   !projName xorig yorig xcell ycell nrows ncols
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
        call system('rm tmp_gdal.txt')
 end subroutine


end program



!ESTRUCTURA GENERAL de MEGAN_PREP_Code_Jan_2022:
!(1) Leen Namelist
!(2) call wrf_file         (leen wrfout para sacar algunos parametros de grilla y proyeccion)
!(3) call  megan2_bioemiss (creo que aca es donde interpola)
!(4) call write_*
!       prepmegan4cmaq_arid.f90:         write(99,'(a)')"CID, ICELL, JCELL,  ARID"
!       prepmegan4cmaq_cantype.f90:      write(99,'(a)')"CID,ICELL,JCELl,NEEDL,TROPI,BROAD,SHRUB,HERB,CROP"
!       prepmegan4cmaq_ecotype.f90:      write(99,'(a)')"gridID,EcotypeID,EcotypeFrac"
!       prepmegan4cmaq_ef.f90:           write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,DSRAD,DTEMP,ISOP,MYRC,SABI,LIMO,A_3CAR,OCIM,BPIN,APIN,OMTP,FARN,BCAR,OSQT,MBO,MEOH,ACTO,CO,NO,BIDER,STRESS,OTHER"
!       prepmegan4cmaq_fert.f90:         write(99,888)"CELL_ID,X,Y,LAT,LONG",fertday
!       prepmegan4cmaq_grwform.bck.f90:  write(99,'(a)')"gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac"
!       prepmegan4cmaq_grwform.f90:      write(99,'(a)')"gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac"
!       prepmegan4cmaq_lai.f90:          write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LAI01,LAI02,LAI03,LAI04,LAI05,LAI06,LAI07,LAI08,LAI09,LAI10, &
!       prepmegan4cmaq_landtype.f90:     write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LANDTYP"
!       
!       prepmegan4cmaq_nitrogen.f90:     write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,NITROGEN01,NITROGEN02,NITROGEN03,NITROGEN04,NITROGEN05,NITROGEN06,NITROGEN07,NITROGEN08,NITROGEN09,NITROGEN10, &
!       prepmegan4cmaq_non_arid.f90:     write(99,'(a)')"CID, ICELL, JCELL,  NON_ARID"
!       prepmegan4cmaq_pft.f90:          write(99,'(a)')"CID, ICELL, JCELL,  NT_EG_TEMP,  NT_DC_BORL,  NT_EG_BORL, &
!       prepmegan4cmaq_w126.f90:         write(99,'(a)')"CID, ICELL, JCELL,  W126(ppm-hours)"
!(5) clean-up (deallocate)
