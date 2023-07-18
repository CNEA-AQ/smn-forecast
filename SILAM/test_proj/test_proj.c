#include <stdio.h>
#include <string.h>
#include <proj.h>
//#include <proj_api.h>

#define d_pi 3.14159265359

int main() {
    int i;    
    double drad_to_deg = 180.0/d_pi;
    double ddeg_to_rad = d_pi/180.0;

    char *lonlat_proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
    char *hirlam_rll_proj4 ="+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0";
    char *lcc_example="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45";

    PJ *P;			//define transformation/projection object
    PJ_COORD c[4], c_out[4];	//define coordinates cin and cout

    //ini proj/transformation
    P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,lonlat_proj4,hirlam_rll_proj4,NULL);
    
    //ini coordinates
    c[0].lpzt.lam=-26.0; //*ddeg_to_rad;
    c[0].lpzt.phi=-35.0; //*ddeg_to_rad;
    c[1].lpzt.lam=-26.0; //*ddeg_to_rad;
    c[1].lpzt.phi= 22.5; //*ddeg_to_rad;
    c[2].lpzt.lam= 40.0; //*ddeg_to_rad;
    c[2].lpzt.phi= 22.5; //*ddeg_to_rad;
    c[3].lpzt.lam= 40.0; //*ddeg_to_rad;
    c[3].lpzt.phi=-35.0; //*ddeg_to_rad;


    //(1) (ll_WGS84) -> Hirlam_rll .
    printf("%s\n","rotated-pole coordinates lon, lat");
    printf("%s\n",hirlam_rll_proj4);
    for (i=0;i<=3;i++){

	c_out[i] = proj_trans(P, PJ_IDENT, c[i]);
	printf("%-3d %20.15f %20.15f\n", i+1, c_out[i].lpzt.lam, c_out[i].lpzt.phi);
    }
    //(2) Hirlam_rll -> (ll_WGS84)
    printf("%s\n","Same in WGS84 coordinates lon, lat");
    printf("%s\n",lonlat_proj4);
    for (i=0;i<=3;i++){
    
    	c[i] = proj_trans(P, PJ_INV, c_out[i]);
        printf("%-3d %20.15f %20.15f\n", i+1, c[i].lpzt.lam, c[i].lpzt.phi);
    }
    //(2)((ll_WGS84) -> Hirlam_rll
    printf("%s\n","And back! lon, lat");
    printf("%s\n",hirlam_rll_proj4);
    for (i=0;i<=3;i++){
     	
     	c_out[i] = proj_trans(P, PJ_FWD, c[i]);
        printf("%-3d %20.15f %20.15f\n", i+1, c_out[i].lpzt.lam, c_out[i].lpzt.phi);                                                                      
    }
    printf("%s","Done");

    // Clean up 
    proj_destroy(P);

    return 0;
}

// subroutine test_proj()
//    implicit none
//
//    integer, parameter :: np = 4
//    real(8), dimension(np) :: lats, lons
//    integer :: i
//
//    character (len=*), parameter :: hirlam_rll_proj4 = &
//           & "+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0"
//    character (len=*), parameter :: lcc_example="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45"
//
//print*,"HOLA ESTOY REVISASNDO PROJ4 NOMA... "
//
//    !! corners of EU meteo in projected lon-lat
//    lons(1:np) = (/-26., -26., 40., 40./)*ddeg_to_rad
//    lats(1:np) = (/-35., 22.5, 22.5, -35./)*ddeg_to_rad
//
//    call msg("rotated-pole coordinates lon, lat")
//    call msg(hirlam_rll_proj4)
//    do i=1,np
//      call msg(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
//    enddo
//
//    call lonalt2proj(hirlam_rll_proj4, lons, lats, np, backwards)
//
//    call msg("Same in WGS84 coordinates lon, lat")
//    call msg(lonlat_proj4)
//    do i=1,np
//      call msg(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
//    enddo
//
//    call lonalt2proj(hirlam_rll_proj4, lons, lats, np, forwards)
//
//    call msg("And back! lon, lat")
//    call msg(hirlam_rll_proj4)
//    do i=1,np
//      call msg(fu_str(i), lons(i)*drad_to_deg,lats(i)*drad_to_deg)
//    enddo
//    call msg("Done")
//
//  end subroutine test_proj
//
//
// HOLISSS...
// HOLA ESTOY REVISASNDO PROJ4 NOMA...
// rotated-pole coordinates lon, lat
// +proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0
// 1  -26.000000000000000       -35.000000000000000
// 2  -26.000000000000000        22.500000000000000
// 3   40.000000000000000        22.500000000000000
// 4   40.000000000000000       -35.000000000000000
// Same in WGS84 coordinates lon, lat
// +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs
// 1  -22.548468112707116        20.537606807248636
// 2  -78.313127936332833        65.570335881774383
// 3   87.834771094127930        53.538491626638915
// 4   33.010219003223881        14.871298797462948
// And back! lon, lat
// +proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0
// 1  -25.999999999999993       -35.000000000000000
// 2  -25.999999999999979        22.500000000000007
// 3   39.999999999999993        22.500000000000000
// 4   39.999999999999993       -35.000000000000000
// Done
//
