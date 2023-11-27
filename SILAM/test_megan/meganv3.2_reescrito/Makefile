.SUFFIXES: .o .f90

FC  =gfortran
LIBS=-L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm
INC =-I/usr/include

FCFLAGS=-O2 -ffree-line-length-none #-Wunused 
LDFLAGS=-O2
OBJS =module_utils.o module_proj.o module_griddesc_reader.o main_megan.o #megvea.o megsea.o
#OBJS =module_utils.o module_proj.o module_griddesc_reader.o module_nc_handler.o module_interp.o megcan.o main_megan.o #megvea.o megsea.o
#OBJS =module_utils.o module_proj.o module_griddesc_reader.o module_nc_handler.o module_interp.o megcan.o megvea.o megsea.o main_megan.o #
EXE  = meganv3.2.exe

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
clean:
	rm -f *.o *.mod                             
