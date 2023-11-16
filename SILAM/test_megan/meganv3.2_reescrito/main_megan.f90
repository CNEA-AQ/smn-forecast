program meganv32
   !
   ! program:     MEGAN
   ! version:     v3.2
   ! description: Biogenic VOCs emission model (by Alex Guenther)
   ! author:      coded by Ramiro Espada
   ! date:        11/2023 
   !
   implicit none
 
   use netcdf   

   use megcan   !canopy calculations
   use megsea   !soil emission activity (SEA): BSNP
   use megvea   !vegetation emision activity (VEA)
   use mgn2mech !species mapping to mechanism

   character(19)  :: start_date, end_date
   character(5)   :: chem_mech='CBM05' !'CBM06' 'CB6A7', 'RACM2','CRACM', 'SAPRC' 'NOCON'
   character(200) :: griddesc_file,gridname,met_file,ct3_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,land_file
   !---
   !Read namelist
   namelist/megan/start_date,end_date,met_file,ct3_file,lai_file,ef_file,ldf_file,ndep_file,fert_file,chem_mech
   read(*,nml=megan, iostat=iostat)
   if( iostat /= 0 ) then
     write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
     stop
   end if
  
   !
   !Initialize stuff



   !do t=1,ntimes

      !
      !get data from meteo, lai, etc.
      call get_data_from_meteo()

      !(este es buen lugar para hacer un broadcast con mpi )

      !
      !calculate canopy values
      call megcan()
 
      !
      !calculate soil (NO) activity
      call megsea() !BDSNP

      !
      !calculate vegetation emission activity
      call megvea()

      !
      !speciate/assign emissions to mechanism 
      call meg2mech()

   !enddo

contains




end program meganv32
