program change_nc_title_value
!
! Programa que le cambia el TITLE a un wrfout por uno que MCIP pueda aceptar.
!
! Para compilarlo:
! > gfortran change_title.f90 -lnetcdf -lnetcdff -I/home/ramiroespada/shared_divqa_data/libs/gcc_6.3.0/netcdf4.4/include -L/home/ramiroespada/shared_divqa_data/libs/gcc_6.3.0/netcdf4.4/lib -o change_title.exe
! 
! Para correrlo:
! > ./change_title.f90 <ruta_a_NetCDF_File_que_deseo_modificar>
!
   use netcdf

  implicit none
  character(300) :: TITLE,oldTITLE,ncfile
  integer :: ncid,var_id,nargs

  nargs = command_argument_count()
  if(nargs/=1) stop "Error, especificar archivo a modificar : > ./programa <ruta_a_NetCDF_File> "
  CALL get_command_argument(1, ncFile)
  TITLE=" OUTPUT FROM WRF V4    MODEL"
  call check(nf90_open(trim(ncfile), nf90_write, ncid ))
       call check(nf90_get_att(ncid, nf90_global, "TITLE", oldTITLE) );
       call check(nf90_put_att(ncid,  nf90_global,"TITLE", trim(TITLE) ))
  call check(nf90_close(ncid))

  print*,"Cambiaste el TITULO al archivo ",trim(ncfile)
  print*," de:", trim(oldTITLE)
  print*,"  a:", trim(TITLE)

contains
subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
      write(*,*) nf90_strerror(status)
      stop 'netcdf error'
   end if
end subroutine check
end program

