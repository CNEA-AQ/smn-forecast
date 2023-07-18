program test_proj

!use, intrinsic :: ISO_C_BINDING, only: c_int,c_char,c_double,c_ptr,c_null_ptr,c_associated,c_size_t,c_loc,c_funptr,c_null_funptr
use ISO_C_BINDING
implicit none

#ifdef USE_PROJ6

type, bind(C) :: PJ_XY
    real(c_double) :: x
    real(c_double) :: y
endtype

interface
  function proj_create_crs_to_crs(s_srs,t_srs) bind(C,name='proj_create_crs_to_crs')
    use iso_c_binding
    implicit none
    type(C_PTR) :: proj_create_crs_to_crs
    character(kind=C_CHAR) :: s_srs,t_srs   !source and target SRS.
  end function

  function proj_trans(pj, dir, coords) bind(C,name='proj_trans')
    use iso_c_binding
    implicit none
    type(C_PTR), value :: PJ
    character(C_CHAR), value :: dir !PJ_FWD, PJ_INV, PJ_IDENT
    type(C_PTR) :: coords
    type(C_PTR) :: proj_trans
  end function

end interface
#endif


!MAIN:
!     character (len=*), parameter :: hirlam_rll_proj4 = &
!           & "+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0"
!    character (len=*), parameter :: lcc_example="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"
!
!
!    !! corners of EU meteo in projected lon-lat
!    lons(1:np) = (/-26., -26., 40., 40./)*ddeg_to_rad
!    lats(1:np) = (/-35., 22.5, 22.5, -35./)*ddeg_to_rad
!
!    print*,("rotated-pole coordinates lon, lat")
!    print*,(hirlam_rll_proj4)
!    do i=1,np
!      print*,(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
!    enddo
!
!    call lonalt2proj(hirlam_rll_proj4, lons, lats, np, backwards)
!
!    print*,("Same in WGS84 coordinates lon, lat")
!    print*,(lonlat_proj4)
!    do i=1,np
!      print*,(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
!    enddo
!
!    call lonalt2proj(hirlam_rll_proj4, lons, lats, np, forwards)
!
!    print*,("And back! lon, lat")
!    print*,(hirlam_rll_proj4)
!    do i=1,np
!      print*,(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
!    enddo
!    print*,("Done")
!

contains

function create_proj(s_srs,t_srs)       result(P)
        implicit none
        type(c_ptr) :: P
        character(len=*) :: s_srs,t_srs
        P=proj_create_crs_to_crs(s_srs,t_srs)
end function

subroutine latlon2proj(P,direction,c,x,y) 
        implicit none
        integer :: i
        type(c_ptr) :: P
        type(PJ_XY) :: c
        type(PJ_XY) :: c_out
        integer :: np
        character(3) :: direction !FWD,INV,ITY
        character(7) :: dir       !PJ_FWD,PJ_INV,PJ_IDENT
        real,intent(inout) :: x,y !output

        if( direction == "FWD") then
                dir="PJ_FWD"
        elseif( direction == "ITY")then
                dir="PJ_IDEN"
        elseif( direction == "BWD")then
                dir="PJ_INV"
        else
                print*, "ERROR: Direction not known";stop;
        endif

        do i=1,1,np
               c_out=proj_trans(P, dir, c_loc(c))
               x=c_out%x
               y=c_out%y
        enddo

end subroutine






end program

