program meganv32
   !
   ! program:     MEGAN  v3.2
   ! description: Biogenic VOCs emission model (by Alex Guenther)
   ! coded:       by Ramiro Espada (date: 11/2023)
    use netcdf   
    use utils_mod         !utils
    use proj_mod          !coordinate transformation functions
    use nc_handler_mod    !functions to deal with netcdf files
    use interpolate_mod   !function for interpolation/regridding
    use readGRIDDESC_mod  !GRIDDESC reader
 
   !#use megcanopy !canopy calculations
   !#use megsea    !soil emission activity (SEA): BSNP
   !#use megvea    !vegetation emision activity (VEA)
   !#use mgn2mech  !species mapping to mechanism

   implicit none

  type(proj_type) :: proj
  type(grid_type) :: grid

   !namelist variables:
   integer           :: iostat
   integer           :: ncid,var_id
   character(len=19) :: start_date, end_date
   character(5)      :: chem_mech='CBM05' !'CBM06' 'CB6A7', 'RACM2','CRACM', 'SAPRC' 'NOCON'
   logical           :: prep_megan=.false.
   character(250)    :: griddesc_file,gridname,met_file,ct3_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file

   !date vars:
   integer      :: end_date_s, current_date_s
   character(4) :: YYYY
   character(3) :: DDD
   character(2) :: MM,DD,HH,current_day,current_month
   
   !input vars:
   integer, allocatable,dimension(:,:)   :: arid,non_arid,landtyp       !(x,y)   from_prep_megan
   real,    allocatable,dimension(:,:,:) :: ct3,ef,ldf                  !(x,y,*) from prep_megan
   real,    allocatable,dimension(:,:)   :: lai,ndep,fert               !(x,y,t) from prep_megan
   real,    allocatable,dimension(:,:,:) :: tmp,par,u10,v10,pre,hum,smo !(x,y,t) from wrfout
   !output vars:
   character(25)    :: outfile,diagfile

   !---
   !Read namelist
   namelist/megan/start_date,end_date,met_file,ct3_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file,chem_mech,griddesc_file,gridname
   read(*,nml=megan, iostat=iostat)
   if( iostat /= 0 ) then
     write(*,*) 'megan: failed to read namelist; error = ',iostat
     stop
   end if

   !---
   print '("Levanto grid y proj parameters.. ")'
   call read_GRIDDESC(griddesc_file,gridname, proj, grid)

   !----
   print '("Preparo datos estáticos.. ")'
   call prep_static_data(grid,arid,non_arid,landtyp,ct3,ef,ldf)

    current_date_s = atoi( date(start_date, "%s") )
        end_date_s = atoi( date(  end_date, "%s") )
   current_day="99" !inicializo con valor absurdo.
   !LOOP OVER TIME
   do while (current_date_s <= end_date_s)

      YYYY=date("@"//itoa(current_date_s), "%Y")   !año
        MM=date("@"//itoa(current_date_s), "%m")   !mes
        DD=date("@"//itoa(current_date_s), "%d")   !dia
       DDD=date("@"//itoa(current_date_s), "%j")   !dia juliano
        HH=date("@"//itoa(current_date_s), "%H")   !hora
      print '("Día: ",A4,"-",A2,"-",A2," (día: ",A3,"), hora: ",2A,".")',YYYY,MM,DD,DDD,HH
      
      if ( current_day /= DD ) then
         current_day=DD
         call prep_dynamic_daily_data(grid,tmp,par,u10,v10,pre,hum,smo,fert)
            if ( current_month /= MM ) then
               current_month=MM
               call prep_dynamic_monthly_data(grid,MM,lai,ndep)
            endif    
      endif    

     !# call megcan(datetime(t), lat, lon,                     &!calculate canopy tmp and rad parameters
     !#        tmp,par,wind,pres,qv,ctf,lai,                   &  !(INPs)
     !#        ShadeleafTK,SunleafTK,SunFrac,SunPPFD,ShadePPFD )  !(OUTs)

     !# call megsea(datetime(t),lat,                           &!calculate soil (NO) activity
     !#         SLTYP,CTF,LAIc,                                &  !Soil type, CT3, LAIc
     !#         TMP,SOILM1,SOILM2,SOILT,PRECADJ,               &  !temp, soil moisture, soil temp, precip adjustment
     !#         CFNO, CFNOG, GAMSM,                            &  !Outs: Emiss activity (EA) of Crop & EA Grass, Soil moisture for isoprene
     !#         GAMNO, BDSNP_NO                                )  !Outs: Final NO emission activity BDSNP NO emissions(nmol/s/m2)

     !# call megvea(datetime(t),                               &!calculate vegetation emission activity
     !#         LAYERS,                                        &  !nlayers
     !#         LAIp, LAIc, LDF_in,                            &  !LAI (past), LAI (current), LDF
     !#         GAMSM_in, MaxT, MinT, MaxWS, D_TEMP, D_PPFD,   &  !
     !#         SunFrac,SUNT, SHAT,SUNP,SHAP,                  &  !from "MEGCAN"
     !#         EMIS_RATE, NON_DIMGARMA                        )  !Outs

     !#!call write_diagnostic_file()
     !#call write_output_file()

     current_date_s=current_date_s + 3600  !siguiente hora!
   enddo

   !#!speciate/assign emissions to mechanism 
   !#call meg2mech()

contains

subroutine prep_static_data(g,arid,non_arid,landtyp,ct3,ef,ldf)
  implicit none
  type(grid_type) :: g
  integer :: k, NRTYP
  integer, allocatable,dimension(:,:)  ,intent(inout) :: arid,non_arid,landtyp       !(x,y)   from_prep_megan
  real,    allocatable,dimension(:,:,:),intent(inout) :: ct3,ef,ldf                  !(x,y,*) from prep_megan

  character(len=10),dimension(19) :: ef_vars=["EF_ISOP   ", "EF_MBO    ", "EF_MT_PINE", "EF_MT_ACYC", "EF_MT_CAMP", "EF_MT_SABI", "EF_MT_AROM", "EF_NO     ", "EF_SQT_HR ", "EF_SQT_LR ", "EF_MEOH   ", "EF_ACTO   ", "EF_ETOH   ", "EF_ACID   ", "EF_LVOC   ", "EF_OXPROD ", "EF_STRESS ", "EF_OTHER  ", "EF_CO     "]
  character(len=5),dimension(4)   :: ldf_vars=["LDF03","LDF04","LDF05","LDF06"]
  NRTYP=6 !(6 canopy types)
  !@ prep static data:
  !Allocation of variables to use
  allocate(    arid(g%nx,g%ny               ))
  allocate(non_arid(g%nx,g%ny               ))
  allocate( landtyp(g%nx,g%ny               ))
  allocate(     ct3(g%nx,g%ny,NRTYP         ))
  allocate(      ef(g%nx,g%ny,size( ef_vars)))
  allocate(     ldf(g%nx,g%ny,size(ldf_vars)))
                                                                                
  LANDTYP =get_2d_var_from_nc_int(land_file,'LANDTYPE',g%nx, g%ny)
  ARID    =get_2d_var_from_nc_int(land_file,'ARID'    ,g%nx, g%ny)
  NON_ARID=get_2d_var_from_nc_int(land_file,'NONARID' ,g%nx, g%ny)
  !!!CT3     =get_3d_var_from_nc(     ct3_file,'CTS'     ,g%nx, g%ny, NRTYP)

  do k=1,size(ef_vars)
     EF(:,:,k)=get_2d_var_from_nc(  ef_file, trim( ef_vars(k)), g%nx, g%ny) 
  enddo
  do k=1,size(ldf_vars)
     LDF(:,:,k)=get_2d_var_from_nc(ldf_file, trim(ldf_vars(k)), g%nx, g%ny) 
  enddo
end subroutine


subroutine prep_dynamic_daily_data(g,tmp,par,u10,v10,pre,hum,smo,fert)
  implicit none
  type(grid_type) :: g
  integer :: ncid,time_dimid,time_len
  real, allocatable,dimension(:,:)    , intent(inout) :: fert
  real, allocatable,dimension(:,:,:), intent(inout) :: u10,v10,tmp,par,pre,hum,smo
  
  !check how many hours is in meteo file
  call check(nf90_open(trim(met_file), nf90_write, ncid     ))
    call check (nf90_inq_dimid(ncid, "Time", time_dimid     ))
    !call check (nf90_inq_dimlen(ncid,time_dimid,time_len))
    call check (nf90_inquire_dimension(ncid,time_dimid,len=time_len))
  call check(nf90_close(ncid))
  print*,"time len: ",time_len
  time_len=time_len-1
  allocate( tmp(g%nx,g%ny,time_len))
  allocate( par(g%nx,g%ny,time_len))
  allocate( u10(g%nx,g%ny,time_len))
  allocate( v10(g%nx,g%ny,time_len))
  allocate( pre(g%nx,g%ny,time_len))
  allocate( hum(g%nx,g%ny,time_len))
  allocate( smo(g%nx,g%ny,time_len))
                                                                                                           
  U10 =get_3d_var_From_nc(met_file,'U10'   , g%nx, g%ny, time_len)
  V10 =get_3d_var_From_nc(met_file,'V10'   , g%nx, g%ny, time_len)
  TMP =get_3d_var_From_nc(met_file,'T2'    , g%nx, g%ny, time_len)
  PAR =get_3d_var_From_nc(met_file,'SWDOWN', g%nx, g%ny, time_len)
  PRE =get_3d_var_From_nc(met_file,'PSFC'  , g%nx, g%ny, time_len)
  HUM =get_3d_var_From_nc(met_file,'Q2'    , g%nx, g%ny, time_len)
  SMO =get_3d_var_From_nc(met_file,'SMOIS' , g%nx, g%ny, time_len)

  allocate( fert(g%nx,g%ny ))
  FERT=get_2d_var_from_nc(fert_file,'FERT'//DDD, g%nx, g%ny)

end subroutine
 
subroutine prep_dynamic_monthly_data(g,MM,lai,ndep)
  implicit none
  type(grid_type) :: g
  character(len=2) :: MM
  real, allocatable,dimension(:,:), intent(inout) :: lai,ndep

  allocate( lai(g%nx,g%ny)) 
  allocate(ndep(g%nx,g%ny))
  LAI=get_2d_var_from_nc(  lai_file, 'LAI'//MM     , g%nx, g%ny)
  NDEP=get_2d_var_from_nc(ndep_file, 'NITROGEN'//MM, g%nx, g%ny)
end subroutine

!!@#subroutine write_output_file()
!!@#  implicit none
!!@#
!!@#  !Creo NetCDF file
!!@#   outfile="out_"//DDD//"_"//MM//"_"//HH//".nc"
!!@#   print*,outfile
!!@#   call createNetCDF(trim(outfile),proj,grid,                     &
!!@#             ['U10','V10','PRE','TMP','PAR','HUM','SMO'],         &
!!@#             spread("FLOAT",1,7), spread("FLOAT",1,7), spread("FLOAT",1,7))
!!@#    !print*,"algun problema acá?"
!!@#   !Abro NetCDF outFile
!!@#   call check(nf90_open(trim(outfile), nf90_write, ncid       ))
!!@#       call check(nf90_inq_varid(ncid,'U10',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,U10(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'V10',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,V10(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'PRE',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,PRE(:,:,1)  )) 
!!@#       call check(nf90_inq_varid(ncid,'TMP',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,TMP(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'PAR',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,PAR(:,:,2)))
!!@#       call check(nf90_inq_varid(ncid,'HUM',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,HUM(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'SMO',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,SMO(:,:,1) ))
!!@#   !Cierro NetCDF outFile
!!@#   call check(nf90_close( ncid ))
!!@#end subroutine
!!@#
!!@#subroutine write_diagnostic_file() !write input and intermediate data so i can check everything is allright
!!@#  implicit none
!!@#
!!@#  !Creo NetCDF file
!!@#   outfile="diag_"//DDD//"_"//MM//"_"//HH//".nc"
!!@#   print*,outfile
!!@#   call createNetCDF(trim(outfile),proj,grid,                     &
!!@#             ['U10','V10','PRE','TMP','PAR','HUM','SMO'],         &
!!@#             spread("FLOAT",1,7), spread("FLOAT",1,7), spread("FLOAT",1,7))
!!@#    !print*,"algun problema acá?"
!!@#   !Abro NetCDF outFile
!!@#   call check(nf90_open(trim(outfile), nf90_write, ncid       ))
!!@#       call check(nf90_inq_varid(ncid,'U10',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,U10(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'V10',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,V10(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'PRE',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,PRE(:,:,1)  )) 
!!@#       call check(nf90_inq_varid(ncid,'TMP',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,TMP(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'PAR',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,PAR(:,:,2)))
!!@#       call check(nf90_inq_varid(ncid,'HUM',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,HUM(:,:,1)   ))
!!@#       call check(nf90_inq_varid(ncid,'SMO',var_id))
!!@#       call check(nf90_put_var(ncid, var_id,SMO(:,:,1) ))
!!@#   !Cierro NetCDF outFile
!!@#   call check(nf90_close( ncid ))
!!@#end subroutine
!!@#
end program meganv32
