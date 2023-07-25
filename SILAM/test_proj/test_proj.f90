program test_proj

use,intrinsic:: ISO_C_BINDING, only: c_int, c_char, c_ptr, c_null_ptr, &
        & c_associated, c_size_t, c_loc, c_funptr, c_null_funptr,c_sizeof
implicit none

interface
  function proj_context_create() bind(C,name='proj_context_create')
    use iso_c_binding
    implicit none
    type(c_ptr) :: proj_context_create
  end function
  subroutine proj_context_destroy(ctx) bind(C,name='proj_context_destroy')
    use iso_c_binding
    implicit none
    type(c_ptr),value :: ctx
  end subroutine
  function proj_create_crs_to_crs(ctx,s_srs,t_srs,area) bind(C,name='proj_create_crs_to_crs') 
    use iso_c_binding                                                                         
    implicit none                                                                             
    type(C_PTR)            :: proj_create_crs_to_crs
    character(kind=C_CHAR) :: s_srs(*)
    character(kind=C_CHAR) :: t_srs(*)
    type(C_PTR), value     :: ctx, area
  end function
  subroutine proj_destroy(PJ) bind(C,name='proj_destroy')
    use iso_c_binding
    implicit none
    type(c_ptr), value :: PJ
  end subroutine
  function proj_trans_generic(pj,dir,x,sx,nx,y,sy,ny,z,sz,nz,t,st,nt) bind(C,name='proj_trans_generic')
    use iso_c_binding
    implicit none
    integer(c_int)       :: proj_trans_generic
    type(c_ptr), value   :: pj
    integer(c_int),value :: dir 
    type(c_ptr),value    :: x,y,z,t
    integer(c_long),value:: sx,sy,sz,st
    integer(c_int),value :: nx,ny,nz,nt
  end function
end interface

!MAIN:
integer :: i
character (len=*), parameter ::     lonlat_proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
character (len=*), parameter :: hirlam_rll_proj4="+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0"
character (len=*), parameter ::      lcc_example="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"

integer, parameter :: np=4
real(8)     :: lats(np),lons(np)

type(c_ptr) :: CTX = c_null_ptr
type(c_ptr) ::  PJ = c_null_ptr

!corners of EU meteo in projected lon-lat
lons(1:np) = (/-26., -26.0, 40.0,  40.0/)
lats(1:np) = (/-35.,  22.5, 22.5, -35.0/)

CTX=proj_context_create() 
PJ=proj_create_crs_to_crs(CTX, lonlat_proj4, hirlam_rll_proj4, C_NULL_PTR)
    
   call latlon2proj(PJ, lons, lats, np, "PJ_IDENT")
   
   print*,("rotated-pole coordinates lon, lat")
   print*,(hirlam_rll_proj4)
   do i=1,np
     print*,i, lons(i),lats(i)
   enddo
   
   call latlon2proj(PJ, lons, lats, np, "PJ_INV")
   
   print*,("Same in WGS84 coordinates lon, lat")
   print*,(lonlat_proj4)
   do i=1,np
     print*,i, lons(i),lats(i)
   enddo
   
   call latlon2proj(PJ, lons, lats, np, "PJ_FWD")
   
   print*,("And back! lon, lat")
   print*,(hirlam_rll_proj4)
   do i=1,np
     print*,i, lons(i),lats(i)
   enddo

   print*,("Done")

!Clean:
call proj_destroy(PJ)
call proj_context_destroy(CTX)

contains

subroutine latlon2proj(PJ,x,y,np,direction) 
  implicit none
  integer             :: i, status
  type(c_ptr)         :: PJ
  integer, intent(in) :: np
  real(8), dimension (:), target, intent(inout) :: x, y
  character(*)        :: direction !PJ_FWD,PJ_INV,PJ_IDENT
  integer :: dir
  integer(8) :: sp

  if      ( trim(direction) == "PJ_FWD"  ) then
     dir=1;
  else if ( trim(direction) == "PJ_IDENT") then                     
     dir=0;
  else if ( trim(direction) == "PJ_INV"  ) then                   
     dir=-1;
  endif

  sp=c_sizeof(x(1))
  status=proj_trans_generic( PJ, dir,              &
                           & c_loc(x(1)),sp, np,   &
                           & c_loc(y(1)),sp, np,   &
                           & c_null_ptr ,sp, 0,    &
                           & c_null_ptr ,sp, 0)  
end subroutine



end program

