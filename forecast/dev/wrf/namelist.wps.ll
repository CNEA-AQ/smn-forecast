&share
 wrf_core = 'ARW',
 max_dom = 1,
 start_date=2024-06-14_00:00:00,
 end_date  =2024-06-15_00:00:00,
 interval_seconds=21600,
 io_form_geogrid = 2,		    !  1:bin; 2:nc; 3:grb1
 active_grid = .true.,
 opt_output_from_geogrid_path ='./',
 !debug_level = 100,
/

&geogrid
 parent_id            = 1             ! parent_id         =   1,  
 parent_grid_ratio    = 1             ! parent_grid_ratio =   1,
 i_parent_start       = 1             ! i_parent_start    =   1, 
 j_parent_start       = 1             ! j_parent_start    =   1,


 e_we                 = 240           ! e_we              = 150,
 e_sn                 = 320           ! e_sn              = 250,
 dx                   = 0.120         ! dx=18000,
 dy                   = 0.120         ! dy=18000,
 map_proj             = 'lat-lon'     ! map_proj='lambert',
 ref_lat              = -38           ! ref_lat=-38.73,
 ref_lon              = -64           ! ref_lon=-62.23,
 !pole_lat            = 52.004       ! truelat1=-33.721,
 !pole_lon            = 0            ! truelat2=-33.721,
 stand_lon            = 64           ! stand_lon=-59.834,

 
 geog_data_res='default',
 geog_data_path = './WPS_GEOG/',
 opt_geogrid_tbl_path = './',
/
&ungrib
 out_format = 'WPS',
 prefix     = 'UNGRIB_FILE',
 pmin = 100.,
/

&metgrid
 fg_name          = 'UNGRIB_FILE'
 io_form_metgrid  = 2, 
 opt_metgrid_tbl_path = './',
 opt_output_from_metgrid_path = './',
 process_only_bdy = 5,
/
