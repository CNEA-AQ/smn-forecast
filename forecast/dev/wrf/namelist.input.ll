&time_control
	run_days=1,
	start_year=2024,
	start_month=06,
	start_day=12,
	start_hour=00,
 start_minute       = 00,
 start_second       = 00,
	end_year=2024,
	end_month=06,
	end_day=13,
	end_hour=00,
 end_minute         = 00,
 end_second         = 00,
	interval_seconds=21600,
 input_from_file    = .true.,
 history_interval   = 60,
 frames_per_outfile = 72,
 restart            = .false.,
 restart_interval   = 0,
 io_form_history    = 2,
 io_form_restart    = 2,
 io_form_input      = 2,
 io_form_boundary   = 2,
 debug_level        = 1,
 auxinput1_inname   = "met_em.d<domain>.<date>",
/

&domains
	time_step=90,
 time_step_fract_num = 0,
 time_step_fract_den = 1,
 max_dom              = 1           ! max_dom= 1,
 e_we                 = 240         ! e_we=175,
 e_sn                 = 320         ! e_sn=300,
 dx                   = 13341.297   ! dx=16000,
 dy                   = 13341.297   ! dy=16000,
 
 e_vert = 45,
 p_top_requested         = 1000, !5000,
 num_metgrid_soil_levels = 4,
	num_metgrid_levels=34,
 use_levels_below_ground = .true.,
 use_surface             = .true.,
 grid_id                 = 1,
 parent_id               = 0,
 i_parent_start          = 1,
 j_parent_start          = 1,
 parent_grid_ratio       = 1,
 parent_time_step_ratio  = 1,
 feedback                = 0,
 smooth_option           = 0,
 !interp_type             = 1,
 !lagrange_order          = 9, 
 !t_extrap_type           = 2,
 force_sfc_in_vinterp    = 0,
 zap_close_levels        = 50,
/

&bdy_control
 spec_bdy_width                      = 5,
 spec_zone                           = 1,
 relax_zone                          = 4,
 specified                           = .true.,
 nested                              = .false.,
/

&namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
/

&dynamics
 w_damping                           = 0,
 diff_opt                            = 1,
 km_opt                              = 4,
 diff_6th_opt                        = 0,
 diff_6th_slopeopt                   = 1,
 diff_6th_thresh                     = 0.10,
 diff_6th_factor                     = 0.12,
 base_temp                           = 290.,
 damp_opt                            = 3,
 zdamp                               = 5000.,
 dampcoef                            = 0.2,
 khdif                               = 0,
 kvdif                               = 0,
 non_hydrostatic                     = .true.,
 moist_adv_opt                       = 1,
 scalar_adv_opt                      = 1,
 h_mom_adv_order                     = 5,
 v_mom_adv_order                     = 3,
 h_sca_adv_order                     = 5,
 v_sca_adv_order                     = 3,
/


&physics
 mp_physics                          = 6,
 bl_pbl_physics                      = 2,
 sf_sfclay_physics                   = 2,
 ra_lw_physics                       = 4,
 ra_sw_physics                       = 4,
 radt                                = 10,
 sf_surface_physics                  = 4,
 bldt                                = 0,
 cu_physics                          = 0,
 cudt                                = 5,
 isfflx                              = 1,
 ifsnow                              = 1,
 icloud                              = 1,
 surface_input_source                = 1,
 num_soil_layers                     = 4,
 sf_urban_physics                    = 0,
 prec_acc_dt                         = 60,
/


!&chem
!	chem_opt=401,
! chemdt                              = 10,
! 
!	chem_in_opt=0,
! 
! io_style_emissions                  = 0,
! 
!	emiss_opt=0,
!	emiss_inpt_opt=0,
!	kemit=1,
! 
! bio_emiss_opt                       = 0,
! bioemdt                             = 0,
! ne_area                             = 0,
!
! biomass_burn_opt                    = 0,
! plumerisefire_frq                   = 30,
! 
! emiss_opt_vol                       = 0,
! emiss_ash_hgt                       = 20000.,
!
!	dust_opt=1,
!	!dust_schme=0,
!
! seas_opt                            = 0,
! dmsemis_opt                         = 0,
!
! gas_bc_opt                          = 0,
! gas_ic_opt                          = 0,
! aer_bc_opt                          = 1,
! aer_ic_opt                          = 1,
!
! gaschem_onoff                       = 0,
!	aerchem_onoff=1,
! wetscav_onoff                       = 0,
! cldchem_onoff                       = 0,
! vertmix_onoff                       = 1,
! 
! gas_drydep_opt                      = 0,
!	aer_drydep_opt=1,
! depo_fact                           = 0.25,
!
! phot_opt                            = 0,
! photdt                              = 0,
! aer_ra_feedback                     = 0,
!	aer_op_opt=1,
!	opt_pars_out=0,
!/
