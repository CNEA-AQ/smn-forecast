.SUFFIXES: .o .f90 .f .F90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm -lproj -fopenmp
INC  = -I/usr/include
FFLAGS = -O2 -ffree-line-length-none -x f95-cpp-input -DUSE_PROJ6 -g #-Wunused 
SCRIP_OBJS=SCRIP_KindsMod.o SCRIP_IOUnitsMod.o SCRIP_CommMod.o SCRIP_ErrorMod.o SCRIP_InitMod.o SCRIP_BroadcastMod.o SCRIP_ConfigMod.o constants.o SCRIP_NetcdfMod.o grids.o remap_vars.o SCRIP_RemapParticleMod.o timers.o remap_conservative.o remap_read.o remap_bilinear.o remap_distance_weight.o remap_mod.o remap_write.o remap_bicubic.o
OBJS = ${SCRIP_OBJS} proj_silam.o mod_SCRIP.o test_remap.o
EXE  = a.out

%.o: %.F90
	${FC} ${FFLAGS} -c ${INC} $<
%.o: %.f  
	${FC} ${FFLAGS} -c ${INC} $<
%.o: %.f90
	${FC} ${FFLAGS} -c ${INC} $<
finn2cmaq: ${OBJS}
	${FC} -o ${EXE} ${OBJS} ${LIBS} 

clean:
	rm -f *.o *.mod

