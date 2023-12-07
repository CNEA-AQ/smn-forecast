program meganv32
   !
   ! program:     MEGAN  v3.2
   ! description: Biogenic VOCs emission model (by Alex Guenther)
   ! coded:       by Ramiro Espada (date: 11/2023)
    use netcdf   
    use utils_mod         !utils
    use proj_mod          !coordinate transformation functions
    use readGRIDDESC_mod  !GRIDDESC reader
 
   use megan_ini   !initialize and set some global parameters based on namelist options
   use meg_can     !canopy calculations
   use meg_sea     !soil emission activity (SEA): BSNP
   use meg_vea     !vegetation emision activity (VEA)
   !#use meg_nox   !NOx emision activity (VEA)
   use mgn2mech    !species mapping to mechanism

   implicit none

   type(proj_type) :: proj
   type(grid_type) :: grid

   !namelist variables:
   integer           :: iostat
   integer           :: ncid,var_id
   character(len=19) :: start_date, end_date
   character(len=16) :: mechanism='CBM05' !'CBM06' 'CB6A7', 'RACM2','CRACM', 'SAPRC' 'NOCON'
   !logical           :: prep_megan=.false.
   character(250)    :: griddesc_file,gridname,met_file,ctf_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file
   character(4)      :: lsm
   !date vars:
   integer      :: end_date_s,current_date_s
   character(4) :: YYYY
   character(3) :: DDD
   character(2) :: MM,DD,HH,current_day,current_month
   
   !input variables:
   real,    allocatable, dimension(:,:)     :: lon,lat                     !(x,y)   from wrfout
   integer, allocatable, dimension(:,:)     :: arid,non_arid,landtype      !(x,y)   from prep_megan
   real,    allocatable, dimension(:,:,:)   :: ef,ldf_in                   !(x,y,*) from prep_megan
   real,    allocatable, dimension(:,:,:,:) :: ctf                         !(x,y,*) from prep_megan
   !real,    allocatable, dimension(:,:,:)  :: ctf                         !(x,y,*) from prep_megan
   !real,    allocatable, dimension(:,:)    :: needl,tropi,broad,shrub,grass,crop
   real,    allocatable, dimension(:,:)     :: lai,ndep,fert               !(x,y,t) from prep_megan
   real,    allocatable, dimension(:,:,:)   :: tmp,rad,u10,v10,pre,hum,smois !(x,y,t) from wrfout
   integer, allocatable, dimension(:,:,:)   :: styp                        !(x,y,t) from wrfout
   !intermediate vars:
   integer :: t!,i,j,k
   integer,parameter :: layers=5
   real,    allocatable, dimension(:,:,:)   :: wind                            !(x,y,t) from wrfout
   real,    allocatable, dimension(:,:,:)   :: sunt,shat,sunfrac,sunp,shap     !(x,y,z) ShadeleafTK,SunleafTK,SunFrac,SunPPFD,ShadePPFD from MEGCAN
   real,    allocatable, dimension(:,:)     :: gamma_sm                        !(x,y) from MEGSEA
   real,    allocatable, dimension(:,:)     :: tmp_min,tmp_max,wind_max,tmp_avg,rad_avg !(x,y) daily meteo vars
   real,    allocatable, dimension(:,:)     :: emis_rate                       !(x,y)        from megvea 
   real,    allocatable, dimension(:,:,:)   :: non_dimgarma                    !(x,y,nclass) from megvea 
     !#   cfno, cfnog,gamma_no, bdsnp_no     !Outs: Final NO emission activity BDSNP NO emissions(nmol/s/m2)
     !#   maxt, mint, maxws, d_temp, d_ppfd, 
   !output vars:
   real,    allocatable, dimension(:,:,:,:) :: out_buffer                      !(x,y,nclass,t)

   !character(25)    :: out_file,diag_file

   !---
   !print '("Leo namelist.. ")'
   namelist/megan/start_date,end_date,met_file,ctf_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file,mechanism,griddesc_file,gridname,lsm
   read(*,nml=megan, iostat=iostat)
   if( iostat /= 0 ) then
     write(*,*) 'megan: failed to read namelist; error = ',iostat
     stop
   end if
   !---
   print '("Levanto grid y proj parameters.. ")'
   call read_GRIDDESC(griddesc_file,gridname, proj, grid)

   !---
   print '("Preparo mecanísmo y especies.. ")'
   call megan_map(mechanism)
   print*,"Mecanísmo: ",mechanism
   print*,"Especies: ",megan_names

   !--- 
   print '("Preparo datos estáticos.. ")'
   call prep_static_data(proj,grid,lat,lon,arid,non_arid,landtype,ctf,ef,ldf_in)
 
   !--- 
   print '("Preparo variables intermedias.. ")' 
   allocate(   shat(grid%nx,grid%ny,layers))!ShadeleafTK
   allocate(   sunt(grid%nx,grid%ny,layers))!ShSunleafTK
   allocate(sunfrac(grid%nx,grid%ny,layers))!  ShSunFrac
   allocate(   sunp(grid%nx,grid%ny,layers))!  ShSunPPFD
   allocate(   shap(grid%nx,grid%ny,layers))!ShShadePPFD
   
   allocate(gamma_sm(grid%nx,grid%ny))       !GAMSM
  
   allocate(   emis_rate(grid%nx,grid%ny))        !for megvea  
   allocate(non_dimgarma(grid%nx,grid%ny,nclass)) !for megvea

   allocate(out_buffer(grid%nx,grid%ny,n_spca_spc,24)) !main out array

   !--- 
   print '("Comienza loop.. ")'
   current_date_s = atoi( date(start_date, "%s"))
       end_date_s = atoi( date(  end_date, "%s"))
   current_day="99" !inicializo con valor absurdo.
   current_month="99"
   do while (current_date_s <= end_date_s)
      YYYY=date("@"//itoa(current_date_s), "%Y")  !año
      MM=date("@"//itoa(current_date_s), "%m")  !mes
      DD=date("@"//itoa(current_date_s), "%d")  !dia
      DDD=date("@"//itoa(current_date_s), "%j")  !dia juliano
      HH=date("@"//itoa(current_date_s), "%H")  !hora
      print '("  Día: ",A4,"-",A2,"-",A2," (día: ",A3,"), hora: ",2A,".")',YYYY,MM,DD,DDD,HH

      t=atoi(HH) !time index
      
      if ( current_day /= DD .or. current_date_s == end_date_s ) then
         if ( current_day /= "99" ) then
            print '("  Escribiendo archivo diario. ")'
            call write_output_file(grid,YYYY,DDD,MECHANISM)
         endif
         current_day=DD
         !call prep_dynamic_daily_data(grid,t,tmp,rad,u10,v10,pre,hum,smois,fert)
         call prep_dynamic_daily_data(grid,tmp,rad,u10,v10,pre,hum,smois,fert)
         if ( current_month /= MM ) then
            current_month=MM
            call prep_dynamic_monthly_data(grid,MM,lai,ndep)
         endif    
      endif    

print*,"==>MEGCAN<=="
      call megcan(atoi(yyyy),atoi(ddd),atoi(hh),                     &!calculate canopy tmp and rad parameters
             grid%nx,grid%ny,layers,                                 & !dimensions
             lat,lon,ctf(:,:,1,:)*0.01,lai,                          & !INPs: 
             tmp(:,:,t),rad(:,:,t),wind(:,:,t),pre(:,:,t),hum(:,:,t),& !INPs: 
             sunt,shat,sunfrac,sunp,shap                             ) !OUTs:  perfiles de temp y flujo fotonico para distintos niveles del canopeo
print*,"==>MEGSEA<=="
      call megsea(                                                   &!calculate soil (NO) activity
             grid%nx,grid%ny,                                        & !dimensions
             lsm,styp(:,:,t),smois(:,:,t),                           & !INPs: land surface model, soil type, soil moisture
             gamma_sm                                                ) !OUTs: soil moisture activity for isoprene

      !soil NO model:
      !call megnox(atoi(yyyy),atoi(ddd),atoi(hh), 
      !       gamma_no, bdsnp_no                                          ) !Outs: Final NO emission activity BDSNP NO emissions(nmol/s/m2)
      !ctf(:,:,1,:),lai,                               & !Soil type, CTF, LAIc
      !tmp(:,:,t),soilm1(:,:,t),soilm2(:,:,t),soilt(:,:,t),precadj(:,:,t),  &  !temp, soil moisture, soil temp, precip adjustment
      !cfno, cfnog, gamma_sm)                                      & !Outs: Emiss activity (EA) of Crop & EA Grass, Soil moisture for isoprene
print*,"==>MEGVEA<=="
      call megvea(                                     &!calculate vegetation emission activity
             grid%nx,grid%ny,layers,                   &  !dimensions
             lai, lai, ldf_in,                         &  !LAI (past), LAI (current), LDF
             gamma_sm,                                 &  !GAMSM: EA response to soil moisture,
             tmp_max,tmp_min,wind_max,tmp_avg,rad_avg, &  !max temp, min temp, max wind, daily temp and daily ppfd
             sunt,shat,sunfrac,sunp,shap,              &  !from "MEGCAN"
             emis_rate,non_dimgarma                    )  !Outs
     !laip=laic
     !----
     !speciate/assign emissions to mechanism 
print*,"==>MGN2MECH<=="
     call convert2mech(grid%nx,grid%ny,  &
                  ef,non_dimgarma,   &
                  out_buffer(:,:,:,t))

     !print '("  Escribiendo archivo diagnóstico.. ")'
     call write_diagnostic_file(proj,grid)

     current_date_s=current_date_s + 3600  !siguiente hora!
   enddo

   print '("Fin. ")'
contains

subroutine prep_static_data(p,g,lat,lon,arid,non_arid,landtype,ctf,ef,ldf_in)
  implicit none
  type(grid_type) :: g
  type(proj_type) :: p
  integer         :: i,j,k,NRTYP
  real,   allocatable,dimension(:,:)    ,intent(inout) :: lat,lon                 !(x,y)   
  integer,allocatable,dimension(:,:)    ,intent(inout) :: arid,non_arid,landtype  !(x,y)   from_prep_megan
  real,   allocatable,dimension(:,:,:)  ,intent(inout) :: ef,ldf_in                  !(x,y,*) from prep_megan
  real,   allocatable,dimension(:,:,:,:),intent(inout) :: ctf                     !(x,y,*) from prep_megan
  !real,   allocatable, dimension(:,:)  :: needl,tropi,broad,shrub,grass,crop

  character(len=10),dimension(19) :: ef_vars=["EF_ISOP   ", "EF_MBO    ", "EF_MT_PINE", "EF_MT_ACYC", "EF_MT_CAMP", "EF_MT_SABI", "EF_MT_AROM", "EF_NO     ", "EF_SQT_HR ", "EF_SQT_LR ", "EF_MEOH   ", "EF_ACTO   ", "EF_ETOH   ", "EF_ACID   ", "EF_LVOC   ", "EF_OXPROD ", "EF_STRESS ", "EF_OTHER  ", "EF_CO     "]
  character(len=5),dimension(4)   :: ldf_vars=["LDF03","LDF04","LDF05","LDF06"]

  NRTYP=6 !(6 canopy types)

  !Allocation of variables to use
  allocate(    arid(g%nx,g%ny               ))
  allocate(non_arid(g%nx,g%ny               ))
  allocate(landtype(g%nx,g%ny               ))
  allocate(      ef(g%nx,g%ny,size( ef_vars)))
  allocate(     ldf_in(g%nx,g%ny,size(ldf_vars)))
  allocate(     ctf(g%nx,g%ny,1,NRTYP       ))
  print*,shape(ctf)
  !LAND                                                                               
  call check(nf90_open(trim(land_file), nf90_write, ncid ))
     call check( nf90_inq_varid(ncid,'LANDTYPE', var_id )); call check( nf90_get_var(ncid, var_id, LANDTYPE ))
     call check( nf90_inq_varid(ncid,'ARID'    , var_id )); call check( nf90_get_var(ncid, var_id, ARID     ))
     call check( nf90_inq_varid(ncid,'NONARID' , var_id )); call check( nf90_get_var(ncid, var_id, NON_ARID ))
  call check(nf90_close(ncid))

  !CTS (canopy type fractions)
  call check(nf90_open(trim(ctf_file), nf90_write, ncid ))
      call check( nf90_inq_varid(ncid,'CTS' , var_id )); call check( nf90_get_var(ncid, var_id, CTF ))
      !call check( nf90_inq_varid(ncid,'NEEDL', var_id )); call check( nf90_get_var(ncid, var_id, NEEDL ))
      !call check( nf90_inq_varid(ncid,'TROPI', var_id )); call check( nf90_get_var(ncid, var_id, TROPI ))
      !call check( nf90_inq_varid(ncid,'BROAD', var_id )); call check( nf90_get_var(ncid, var_id, BROAD ))
      !call check( nf90_inq_varid(ncid,'SHRUB', var_id )); call check( nf90_get_var(ncid, var_id, SHRUB ))
      !call check( nf90_inq_varid(ncid,'GRASS', var_id )); call check( nf90_get_var(ncid, var_id, GRASS ))
      !call check( nf90_inq_varid(ncid,'CROP' , var_id )); call check( nf90_get_var(ncid, var_id, CROP  ))
  call check(nf90_close(ncid))
  
  !EF:
  call check(nf90_open(trim(  ef_file), nf90_write, ncid ))
  do k=1,size(ef_vars)
     call check( nf90_inq_varid(ncid,trim( ef_vars(k)), var_id )); call check( nf90_get_var(ncid, var_id , EF(:,:,k) ))
  enddo
  call check(nf90_close(ncid))

  !LDF:
  call check(nf90_open(trim( ldf_file), nf90_write, ncid ))
  do k=1,size(ldf_vars)
     call check( nf90_inq_varid(ncid,trim(ldf_vars(k)), var_id )); call check( nf90_get_var(ncid, var_id ,ldf_in(:,:,k) ))
  enddo
  call check(nf90_close(ncid))

  !LAT & LON
  allocate( lat(g%nx,g%ny ))
  allocate( lon(g%nx,g%ny ))
  do j=1,g%ny
     do i=1,g%nx
        call xy2ll(p,g%xmin + g%dx*i, g%ymin + g%dx*j,lon(i,j),lat(i,j))
     enddo
  enddo
 
end subroutine


!subroutine prep_dynamic_daily_data(g,ini_h,tmp,rad,u10,v10,pre,hum,SMOIS,fert)
subroutine prep_dynamic_daily_data(g,tmp,rad,u10,v10,pre,hum,SMOIS,fert)
  implicit none
  type(grid_type) :: g
  !integer,intent(in) :: ini_h
  integer            :: ini_h,end_h
  integer :: ncid,time_dimid,time_len
  real, allocatable,dimension(:,:)    , intent(inout) :: fert
  real, allocatable,dimension(:,:,:), intent(inout) :: u10,v10,tmp,rad,pre,hum,smois
  character(len=19), dimension(:), allocatable :: Times

print*,"  Preparo inputs dinámicos diarios.."
  !2d arrays:
  if (.not. allocated(fert)) then; allocate( fert(g%nx,g%ny));endif
  call check(nf90_open(trim(fert_file), nf90_write, ncid ))
        call check(   nf90_inq_varid(ncid,'FERT'//DDD, var_id )); call check( nf90_get_var(ncid, var_id , FERT ))
  call check(nf90_close(ncid))
 
  !3d arrays:
  !check how many hours is in meteo file
  call check(nf90_open(trim(met_file), nf90_write, ncid     ))
    call check (nf90_inq_dimid(ncid, "Time", time_dimid     ))
    call check (nf90_inquire_dimension(ncid,time_dimid,len=time_len))
    allocate(Times(time_len))
    call check( nf90_inq_varid(ncid,'Times' , var_id)); call check(nf90_get_var(ncid, var_id,  Times(:) ))
print*,Times
print*,Times(time_len-1)
    read(Times(1         )(12:13), *) ini_h
    read(Times(time_len-1)(12:13), *) end_h

  call check(nf90_close(ncid))
  deallocate(Times)

  !end_h=min(modulo(24-time_len,24), 24)
print*,"ini_h,end_h,time_len: ",ini_h,end_h,time_len
  if (.not. allocated(tmp)  ) then; allocate(  tmp(g%nx,g%ny,24));endif 
  if (.not. allocated(rad)  ) then; allocate(  rad(g%nx,g%ny,24));endif
  if (.not. allocated(u10)  ) then; allocate(  u10(g%nx,g%ny,24));endif
  if (.not. allocated(v10)  ) then; allocate(  v10(g%nx,g%ny,24));endif
  if (.not. allocated(pre)  ) then; allocate(  pre(g%nx,g%ny,24));endif
  if (.not. allocated(hum)  ) then; allocate(  hum(g%nx,g%ny,24));endif
  if (.not. allocated(SMOIS)) then; allocate(smois(g%nx,g%ny,24));endif
  if (.not. allocated(styp) ) then; allocate( styp(g%nx,g%ny,24));endif
  call check(nf90_open(trim(met_file), nf90_write, ncid ))
     call check( nf90_inq_varid(ncid,'U10'   , var_id)); call check(nf90_get_var(ncid, var_id,  U10(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'V10'   , var_id)); call check(nf90_get_var(ncid, var_id,  V10(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'T2'    , var_id)); call check(nf90_get_var(ncid, var_id,  TMP(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'SWDOWN', var_id)); call check(nf90_get_var(ncid, var_id,  RAD(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'PSFC'  , var_id)); call check(nf90_get_var(ncid, var_id,  PRE(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'Q2'    , var_id)); call check(nf90_get_var(ncid, var_id,  HUM(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'SMOIS' , var_id)); call check(nf90_get_var(ncid, var_id,SMOIS(:,:,ini_h:end_h) ))
     call check( nf90_inq_varid(ncid,'ISLTYP', var_id)); call check(nf90_get_var(ncid, var_id, STYP(:,:,ini_h:end_h) ))
  call check(nf90_close(ncid))

  if (.not. allocated(wind)) then; allocate(wind(g%nx,g%ny,24) );endif 
  WIND=SQRT(U10*U10 + V10*V10)
  
  !Ground Incident Radiation [W m-2] to PPFD (Photosynthetic Photon Flux Density [W m-2])
  rad=rad*4.5*0.45 ! ppfd = par   * 4.5     !par to ppfd
                   ! par  = rgrnd * 0.45   !total rad to Photosyntetic Active Radiation (PAR)
  !Dayly variables:
  if (.not. allocated(rad_avg) ) then; allocate( rad_avg(g%nx,g%ny) );endif 
  if (.not. allocated(tmp_avg) ) then; allocate( tmp_avg(g%nx,g%ny) );endif 
  if (.not. allocated(tmp_min )) then; allocate( tmp_min(g%nx,g%ny) );endif 
  if (.not. allocated(tmp_max )) then; allocate( tmp_max(g%nx,g%ny) );endif 
  if (.not. allocated(wind_max)) then; allocate(wind_max(g%nx,g%ny) );endif 
  tmp_min = minval(tmp, dim=3,mask=tmp>0)
  tmp_max = maxval(tmp, dim=3)
  wind_max= maxval(wind,dim=3)
  tmp_avg =   sum(tmp, dim=3)/time_len
  rad_avg =   sum(rad, dim=3)/time_len
end subroutine
 
subroutine prep_dynamic_monthly_data(g,MM,lai,ndep)
  implicit none
  type(grid_type)   :: g
  character(len=2)  :: MM
  real, allocatable,dimension(:,:), intent(inout) :: lai,ndep

print*,"  Preparo inputs dinámicos mensuales.."
  if (.not. allocated(lai) ) then; allocate( lai(g%nx,g%ny));endif
  if (.not. allocated(ndep)) then; allocate(ndep(g%nx,g%ny));endif
  call check(nf90_open( trim(lai_file), nf90_write, ncid ))
     call check( nf90_inq_varid(ncid,'LAI'//MM    , var_id )); call check( nf90_get_var(ncid, var_id , LAI  ))
  call check(nf90_close(ncid))
  call check(nf90_open( trim(ndep_file), nf90_write, ncid ))
     call check( nf90_inq_varid(ncid,'NITROGEN'//MM, var_id )); call check( nf90_get_var(ncid, var_id , NDEP ))
  call check(nf90_close(ncid))
end subroutine


subroutine write_diagnostic_file(p,g) !write input and intermediate data so i can check everything is allright
  implicit none
  type(grid_type) , intent(in) :: g
  type(proj_type) , intent(in) :: p
  character(len=20) :: diag_file
  integer           :: ncid,t_dim_id,x_dim_id,y_dim_id,z_dim_id,s_dim_id,var_dim_id
  integer           :: k!,i,j
  integer           :: ntimes
  character(len=5),dimension(31) :: var_names, var_types, var_dimen

  ntimes=size(U10(1,1,:))

  !Creo NetCDF file
  diag_file="diag_"//YYYY//"_"//DDD//"_"//HH//".nc"

  print*,'   Escrbibiendo diag_file: ',diag_file

data var_names( 1),var_types( 1),var_dimen( 1) / 'LAT  ', 'FLOAT', 'XY   '/
data var_names( 2),var_types( 2),var_dimen( 2) / 'LON  ', 'FLOAT', 'XY   '/
data var_names( 3),var_types( 3),var_dimen( 3) / 'EF   ', 'FLOAT', 'XY   '/
data var_names( 4),var_types( 4),var_dimen( 4) / 'LDF  ', 'FLOAT', 'XY   '/
data var_names( 5),var_types( 5),var_dimen( 5) / 'LAI  ', 'FLOAT', 'XY   '/
data var_names( 6),var_types( 6),var_dimen( 6) / 'NDEP ', 'FLOAT', 'XY   '/
data var_names( 7),var_types( 7),var_dimen( 7) / 'FERT ', 'FLOAT', 'XY   '/
data var_names( 8),var_types( 8),var_dimen( 8) / 'ARID ', 'INT  ', 'XY   '/
data var_names( 9),var_types( 9),var_dimen( 9) / 'NARID', 'INT  ', 'XY   '/
data var_names(10),var_types(10),var_dimen(10) / 'LTYPE', 'INT  ', 'XY   '/
data var_names(11),var_types(11),var_dimen(11) / 'U10  ', 'FLOAT', 'XYT  '/
data var_names(12),var_types(12),var_dimen(12) / 'V10  ', 'FLOAT', 'XYT  '/
data var_names(13),var_types(13),var_dimen(13) / 'PRE  ', 'FLOAT', 'XYT  '/
data var_names(14),var_types(14),var_dimen(14) / 'TMP  ', 'FLOAT', 'XYT  '/
data var_names(15),var_types(15),var_dimen(15) / 'RAD  ', 'FLOAT', 'XYT  '/
data var_names(16),var_types(16),var_dimen(16) / 'HUM  ', 'FLOAT', 'XYT  '/
data var_names(17),var_types(17),var_dimen(17) / 'SMOIS', 'FLOAT', 'XYT  '/
data var_names(18),var_types(18),var_dimen(18) / 'ISTYP', 'INT  ', 'XYT  '/
data var_names(19),var_types(19),var_dimen(19) / 'T_SHA', 'FLOAT', 'XYZ  '/
data var_names(20),var_types(20),var_dimen(20) / 'T_SUN', 'FLOAT', 'XYZ  '/
data var_names(21),var_types(21),var_dimen(21) / 'SunFr', 'FLOAT', 'XYZ  '/
data var_names(22),var_types(22),var_dimen(22) / 'R_SHA', 'FLOAT', 'XYZ  '/
data var_names(23),var_types(23),var_dimen(23) / 'R_SUN', 'FLOAT', 'XYZ  '/
data var_names(24),var_types(24),var_dimen(24) / 'GAMSM', 'FLOAT', 'XY   '/
data var_names(25),var_types(25),var_dimen(25) / 'T_MIN', 'FLOAT', 'XY   '/
data var_names(26),var_types(26),var_dimen(26) / 'T_MAX', 'FLOAT', 'XY   '/
data var_names(27),var_types(27),var_dimen(27) / 'W_MAX', 'FLOAT', 'XY   '/
data var_names(28),var_types(28),var_dimen(28) / 'T_AVG', 'FLOAT', 'XY   '/
data var_names(29),var_types(29),var_dimen(29) / 'R_AVG', 'FLOAT', 'XY   '/
data var_names(30),var_types(30),var_dimen(30) / 'ER   ', 'FLOAT', 'XY   '/
data var_names(31),var_types(31),var_dimen(31) / 'NONDI', 'FLOAT', 'XYC  '/


  call check(nf90_create(trim(diag_file), NF90_CLOBBER, ncid))
    !Defino dimensiones
    call check(nf90_def_dim(ncid, "Time", 1      , t_dim_id       ))
    call check(nf90_def_dim(ncid, "x"   , g%nx   , x_dim_id       ))
    call check(nf90_def_dim(ncid, "y"   , g%ny   , y_dim_id       ))
    call check(nf90_def_dim(ncid, "z"   , layers , z_dim_id       ))
    call check(nf90_def_dim(ncid, "c"   , nclass , s_dim_id       ))
    !Defino variables
    call check(nf90_def_var(ncid,"Times",  NF90_INT    , [t_dim_id], var_id))
    call check(nf90_put_att(ncid, var_id, "units"      , "seconds from file start date" ))
    call check(nf90_put_att(ncid, var_id, "var_desc"   , "date-time variable"      ))
    !Creo variables:
    do k=1, size(var_names)
      if ( trim(var_types(k)) == "FLOAT" ) then
         if ( trim(var_dimen(k)) == 'XY' .or. trim(var_dimen(k)) == 'XYT') then
              call check( nf90_def_var(ncid, trim(var_names(k)) , NF90_FLOAT, [x_dim_id,y_dim_id], var_id)          )
         else if ( trim(var_dimen(k)) == 'XYZ' ) then
              call check( nf90_def_var(ncid, trim(var_names(k)) , NF90_FLOAT, [x_dim_id,y_dim_id,z_dim_id], var_id) )
         else if ( trim(var_dimen(k)) == 'XYC' ) then
              call check( nf90_def_var(ncid, trim(var_names(k)) , NF90_FLOAT, [x_dim_id,y_dim_id,s_dim_id], var_id) )
         endif
      else if ( trim(var_types(k)) == "INT") then   
              if ( trim(var_dimen(k)) == 'XY' .or. trim(var_dimen(k)) == 'XYT'  ) then
              call check( nf90_def_var(ncid, trim(var_names(k)) , NF90_INT  , [x_dim_id,y_dim_id], var_id)          )
         endif
       endif
    end do
  call check(nf90_enddef(ncid))   !End NetCDF define mode

  call check(nf90_open(trim(diag_file), nf90_write, ncid ))
      call check(nf90_inq_varid(ncid,'LAT'  ,var_id )); call check(nf90_put_var(ncid, var_id, LAT        ))
      call check(nf90_inq_varid(ncid,'LON'  ,var_id )); call check(nf90_put_var(ncid, var_id, LON        ))
      call check(nf90_inq_varid(ncid,'EF'   ,var_id )); call check(nf90_put_var(ncid, var_id, EF         ))
      call check(nf90_inq_varid(ncid,'LDF'  ,var_id )); call check(nf90_put_var(ncid, var_id, ldf_in     ))
      !call check(nf90_inq_varid(ncid,'CTS'  ,var_id )); call check(nf90_put_var(ncid, var_id, CTS        ))
      call check(nf90_inq_varid(ncid,'LAI'  ,var_id )); call check(nf90_put_var(ncid, var_id, LAI        ))
      call check(nf90_inq_varid(ncid,'NDEP' ,var_id )); call check(nf90_put_var(ncid, var_id, NDEP       ))
      call check(nf90_inq_varid(ncid,'FERT' ,var_id )); call check(nf90_put_var(ncid, var_id, FERT       ))
      call check(nf90_inq_varid(ncid,'ARID' ,var_id )); call check(nf90_put_var(ncid, var_id, ARID       ))
      call check(nf90_inq_varid(ncid,'NARID',var_id )); call check(nf90_put_var(ncid, var_id, NON_ARID   ))
      call check(nf90_inq_varid(ncid,'LTYPE',var_id )); call check(nf90_put_var(ncid, var_id, LANDTYPE   ))
      call check(nf90_inq_varid(ncid,'U10'  ,var_id )); call check(nf90_put_var(ncid, var_id, U10(:,:,t) ))
      call check(nf90_inq_varid(ncid,'V10'  ,var_id )); call check(nf90_put_var(ncid, var_id, V10(:,:,t) ))
      call check(nf90_inq_varid(ncid,'PRE'  ,var_id )); call check(nf90_put_var(ncid, var_id, PRE(:,:,t) )) 
      call check(nf90_inq_varid(ncid,'TMP'  ,var_id )); call check(nf90_put_var(ncid, var_id, TMP(:,:,t) ))
      call check(nf90_inq_varid(ncid,'RAD'  ,var_id )); call check(nf90_put_var(ncid, var_id, RAD(:,:,t) ))
      call check(nf90_inq_varid(ncid,'HUM'  ,var_id )); call check(nf90_put_var(ncid, var_id, HUM(:,:,t) ))
      call check(nf90_inq_varid(ncid,'SMOIS',var_id )); call check(nf90_put_var(ncid, var_id, SMOIS(:,:,t) ))
      call check(nf90_inq_varid(ncid,'ISTYP',var_id )); call check(nf90_put_var(ncid, var_id,STYP(:,:,t) ))
      call check(nf90_inq_varid(ncid,'T_SHA',var_id )); call check(nf90_put_var(ncid, var_id, SHAT       ))
      call check(nf90_inq_varid(ncid,'T_SUN',var_id )); call check(nf90_put_var(ncid, var_id, SUNT       ))
      call check(nf90_inq_varid(ncid,'SunFr',var_id )); call check(nf90_put_var(ncid, var_id, SunFrac    ))
      call check(nf90_inq_varid(ncid,'R_SHA',var_id )); call check(nf90_put_var(ncid, var_id, SHAP       ))
      call check(nf90_inq_varid(ncid,'R_SUN',var_id )); call check(nf90_put_var(ncid, var_id, SUNP       ))
      call check(nf90_inq_varid(ncid,'GAMSM',var_id )); call check(nf90_put_var(ncid, var_id, GAMMA_SM   ))
      call check(nf90_inq_varid(ncid,'T_MIN',var_id )); call check(nf90_put_var(ncid, var_id, tmp_min   ))
      call check(nf90_inq_varid(ncid,'T_MAX',var_id )); call check(nf90_put_var(ncid, var_id, tmp_max   ))
      call check(nf90_inq_varid(ncid,'T_AVG',var_id )); call check(nf90_put_var(ncid, var_id, tmp_avg   ))
      call check(nf90_inq_varid(ncid,'W_MAX',var_id )); call check(nf90_put_var(ncid, var_id,wind_max   ))
      call check(nf90_inq_varid(ncid,'R_AVG',var_id )); call check(nf90_put_var(ncid, var_id, rad_avg   ))
      call check(nf90_inq_varid(ncid,'ER'   ,var_id )); call check(nf90_put_var(ncid, var_id,emis_rate  ))
      call check(nf90_inq_varid(ncid,'NONDI',var_id )); call check(nf90_put_var(ncid, var_id,non_dimgarma ))
  call check(nf90_close( ncid ))
end subroutine


subroutine write_output_file(g,YYYY,DDD,MECHANISM)
  implicit none
  type(grid_type) , intent(in) :: g
  character(len=5), intent(in) :: MECHANISM
  character(len=4), intent(in) :: YYYY
  character(len=3), intent(in) :: DDD
  !local vars
  integer           :: ncid,t_dim_id,x_dim_id,y_dim_id,z_dim_id,s_dim_id,var_dim_id
  integer           :: k!,i,j
  character(len=25) :: out_file
  !Creo NetCDF file
   out_file="out_"//YYYY//"_"//DDD//"_"//trim(MECHANISM)//".nc"
   print*,"ESCRIBIENDO: ",out_file

   !Crear NetCDF
   call check(nf90_create(trim(out_file), NF90_CLOBBER, ncid))
     !Defino dimensiones
     call check(nf90_def_dim(ncid, "Time", 24     , t_dim_id       ))
     call check(nf90_def_dim(ncid, "x"   , g%nx   , x_dim_id       ))
     call check(nf90_def_dim(ncid, "y"   , g%ny   , y_dim_id       ))
     !Defino variables
     call check(nf90_def_var(ncid,"Times",  NF90_INT    , [t_dim_id], var_id))
     call check(nf90_put_att(ncid, var_id, "units"      , "seconds from file start date" ))
     call check(nf90_put_att(ncid, var_id, "var_desc"   , "date-time variable"      ))
print*,"shape buffer: ",shape(out_buffer)
     !Creo variables:
     do k=1,NMGNSPC
       print*,"Especie: ",k,n_scon_spc,trim(mech_spc(k))
        call check( nf90_def_var(ncid, trim(mech_spc(k)) , NF90_FLOAT, [x_dim_id,y_dim_id,t_dim_id], var_id)   )
     end do
   call check(nf90_enddef(ncid))   !End NetCDF define mode

   !Abro NetCDF y guardo variables de salida
   call check(nf90_open(trim(out_file), nf90_write, ncid       ))
     do k=1,NMGNSPC
       print*,"Especie:2",trim(mech_spc(k))
       call check(nf90_inq_varid(ncid,trim(mech_spc(k)),var_id)); call check(nf90_put_var(ncid, var_id, out_buffer(:,:,k,:) ))
     enddo
   call check(nf90_close( ncid ))
end subroutine write_output_file



end program meganv32

