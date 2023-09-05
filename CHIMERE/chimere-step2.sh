#!/bin/bash
unset LANG

# CHIMERE run in parallel

echo
echo -e "\033[1;43m  Starting run with CHIMERE    \033[0m"
echo "   Model will run as ${nproc_chimere} processes"

#----------------------------------------------
# offline case

if [ ${online} -ne "1" ] && [ ${runwrfonly} == "0" ] ; then

   echo "   offline mode"
   cd ${chimere_tmp}
   if [ -x chimere.e ]; then
      if [ "${machine_name:0:5}" != "irene" ] ; then
         # I Change this to run chimere on neurus and accountneurus and partitionneurus needed to be added to chimere.par
         time srun --verbose -A ${accountneurus} -p ${partitionneurus} -n ${nproc_chimere} ${mpiparams} ./chimere.e
         [ $? -eq 0 ] || { echo "Abnormal termination of chimere.e"; exit 1; }
      else
         time ${my_mpirun} -n ${nproc_chimere} ./chimere.e
         [ $? -eq 0 ] || { echo "Abnormal termination of chimere.e"; exit 1; }
      fi
   else
      echo "No such file ${chimere_tmp}/chimere.e ! Bye."
      exit 1
   fi

fi

#----------------------------------------------
# online case

if [ ${online} -eq "1" ] ; then

   echo "   online mode"

      if [ "${machine_name:0:5}" != "irene" ] ; then
         echo "   running on ${machine_name} with command line :"
         echo "   " time ${my_mpirun} ${mpiparams} -np ${nproc_wrf} ${wrf_exe} : \
                                       -np ${nproc_chimere} ${mpiparams} ./chimere.e
         time ${my_mpirun} ${mpiparams} -np ${nproc_wrf} ${wrf_exe} : \
                                       -np ${nproc_chimere} ${mpiparams} ./chimere.e
         [ $? -eq 0 ] || { echo "Abnormal termination of chimere/WRF"; exit 1; }
      else
         echo "   running on ${machine_name} with ${my_mpirun} :"
         touch chimwrf.conf
         rm chimwrf.conf
         echo ${nproc_wrf} ${wrf_exe} > chimwrf.conf
         echo ${nproc_chimere} ./chimere.e >> chimwrf.conf
         
   echo "Clean up WPS files"
   rm -f ${dirowps}/met_em* met_em* wrfinput_d0? wrfbdy_d0?
   if [ ${do_clean} ==  "full"  ] ; then
      rm -f ${dirowps}/geo_em.*nc
   fi
fi

