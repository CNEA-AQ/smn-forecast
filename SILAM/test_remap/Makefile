.SUFFIXES: .o .f90 .f .F90

FC=gfortran
SE=serial
#SE=mpi
LIBS=-L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm -lproj -fopenmp
INC=-I/usr/include
FFLAGS=-O3 -ffree-line-length-none -x f95-cpp-input -DUSE_PROJ6 #-g -Wunused 

SCRIPDIR=SCRIPT
SCRIP_OBJS=SCRIP_KindsMod.o SCRIP_IOUnitsMod.o SCRIP_ErrorMod.o SCRIP_CommMod.o SCRIP_BroadcastMod.o SCRIP_InitMod.o SCRIP_ConfigMod.o constants.o SCRIP_NetcdfMod.o grids.o remap_vars.o timers.o remap_conservative.o remap_read.o remap_bilinear.o remap_distance_weight.o remap_mod.o remap_write.o remap_bicubic.o # SCRIP_RemapParticleMod.o
OBJS = ${SCRIP_OBJS} proj_silam.o SCRIP_silam.o test_scrip.o
EXE  = a.out

#all:
#	find -maxdepth 1 -type l -delete
#	ls $(SCRIPDIR)/*.F90
#	ls $(SCRIPDIR)/$(SE)/*.F90
#	ln -s $(SCRIPDIR)/$(SE)/*.F90 .
#	ln -s $(SCRIPDIR)/*.f .
#	ln -s $(SCRIPDIR)/*.F90 .
%.o: %.F90
	${FC} ${FFLAGS} -c ${INC} $<
%.o: %.f
	${FC} ${FFLAGS} -c ${INC} $<
%.o: %.f90
	${FC} ${FFLAGS} -c ${INC} $<
$(EXE): ${OBJS}
	${FC} -o ${EXE} $^ ${LIBS} 
clean:
	rm -f *.o *.mod 
clean-links:
	find -maxdepth 1 -type l -delete
