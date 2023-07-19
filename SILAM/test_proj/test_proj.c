#include <stdio.h>
#include <string.h>
#include <proj.h>
//#include <proj_api.h>

#define d_pi 3.14159265359

int main() {
    int i,status;    
    double drad_to_deg = 180.0/d_pi;
    double ddeg_to_rad = d_pi/180.0;

    char *lonlat_proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";
    char *hirlam_rll_proj4 ="+proj=ob_tran +o_proj=longlat +o_lon_p=0 +o_lat_p=30 +lon_0=0";
    char *lcc_example="+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45";

    PJ *P;			//define transformation/projection object
    //PJ_COORD c[4], c_out[4];	//define coordinates cin and cout

    const int np=4;
    double lons[4] = {-26.0,-26.0, 40.0,  40.0};
    double lats[4] = {-35.0, 22.5, 22.5, -35.0};
    //CTX = proj_context_create();	//init context.
    //ini proj/transformation
    //P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,lonlat_proj4,hirlam_rll_proj4,NULL);
    //P = proj_create_crs_to_crs(0,lonlat_proj4,hirlam_rll_proj4,0);
    //P = proj_create_crs_to_crs(0,lonlat_proj4,hirlam_rll_proj4,0);
    P = proj_create_crs_to_crs(NULL,lonlat_proj4,hirlam_rll_proj4,NULL);
    
    ////ini coordinates
    //c[0].lpzt.lam=-26.0; //*ddeg_to_rad;
    //c[0].lpzt.phi=-35.0; //*ddeg_to_rad;
    //c[1].lpzt.lam=-26.0; //*ddeg_to_rad;
    //c[1].lpzt.phi= 22.5; //*ddeg_to_rad;
    //c[2].lpzt.lam= 40.0; //*ddeg_to_rad;
    //c[2].lpzt.phi= 22.5; //*ddeg_to_rad;
    //c[3].lpzt.lam= 40.0; //*ddeg_to_rad;
    //c[3].lpzt.phi=-35.0; //*ddeg_to_rad;
    //
    ////(1) (ll_WGS84) -> Hirlam_rll .
    //printf("%s\n","rotated-pole coordinates lon, lat");
    //printf("%s\n",hirlam_rll_proj4);
    //for (i=0;i<=3;i++){
    //    c_out[i] = proj_trans(P, PJ_IDENT, c[i]);
    //    printf("%-3d %20.15f %20.15f\n", i+1, c_out[i].lpzt.lam, c_out[i].lpzt.phi);
    //}
    ////(2) Hirlam_rll -> (ll_WGS84)
    //printf("%s\n","Same in WGS84 coordinates lon, lat");
    //printf("%s\n",lonlat_proj4);
    //for (i=0;i<=3;i++){
    //
    //	c[i] = proj_trans(P, PJ_INV, c_out[i]);
    //    printf("%-3d %20.15f %20.15f\n", i+1, c[i].lpzt.lam, c[i].lpzt.phi);
    //}
    ////(3)((ll_WGS84) -> Hirlam_rll
    //printf("%s\n","And back! lon, lat");
    //printf("%s\n",hirlam_rll_proj4);
    //for (i=0;i<=3;i++){
    // 	
    // 	c_out[i] = proj_trans(P, PJ_FWD, c[i]);
    //    printf("%-3d %20.15f %20.15f\n", i+1, c_out[i].lpzt.lam, c_out[i].lpzt.phi);                                                                      
    //}
    //printf("%s\n","Done");

    
    printf("================\n%s\n================\n"," WORKARROUND");
    printf("size of double: %ld",sizeof(double));
    //(1) (ll_WGS84) -> Hirlam_rll .
    printf("%s\n","rotated-pole coordinates lon, lat");
    printf("%s\n",hirlam_rll_proj4);
    status=proj_trans_generic ( P, PJ_IDENT, lons, sizeof(double), np, lats, sizeof(double), np, 0, 0, 0, 0, 0, 0 );
    for (i=0;i<np;i++){
        printf("%-3d %20.15f %20.15f\n", i+1, lons[i],lats[i]);                                                                      
    }
    //(2) Hirlam_rll -> (ll_WGS84)
    printf("%s\n","Same in WGS84 coordinates lon, lat");
    printf("%s\n",lonlat_proj4);

    status=proj_trans_generic ( P, -1    , lons, sizeof(double), np, lats, sizeof(double), np, 0, 0, 0, 0, 0, 0 );
    //status=proj_trans_generic ( P, PJ_INV, lons, sizeof(double), np, lats, sizeof(double), np, 0, 0, 0, 0, 0, 0 );
    for (i=0;i<np;i++){
        printf("%-3d %20.15f %20.15f\n", i+1, lons[i],lats[i]);                                                                      
    }
    //(3)((ll_WGS84) -> Hirlam_rll
    printf("%s\n","And back! lon, lat");
    printf("%s\n",hirlam_rll_proj4);

    status=proj_trans_generic ( P,  1    , lons, sizeof(double), np, lats, sizeof(double), np, 0, 0, 0, 0, 0, 0 );
    //status=proj_trans_generic ( P, PJ_FWD, lons, sizeof(double), np, lats, sizeof(double), np, 0, 0, 0, 0, 0, 0 );
    for (i=0;i<np;i++){
        printf("%-3d %20.15f %20.15f\n", i+1, lons[i],lats[i]);                                                                      
    }
    printf("%s\n","Done");

   // Clean up 
   proj_destroy(P);

    return 0;
}
