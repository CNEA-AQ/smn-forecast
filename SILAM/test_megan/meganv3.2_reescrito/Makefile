.SUFFIXES: .o .f90

FC  =gfortran
LIBS=-L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm
INC =-I/usr/include

FCFLAGS=-O2 -ffree-line-length-none #-Wall#-Wunused 
LDFLAGS=-O2

OBJS=module_utils.o module_proj.o module_griddesc_reader.o mgnmap.o megcan.o megvea.o main_megan.o  #megnox.o megsea.o mgn2mech.o 
EXE =meganv3.2.exe

$(EXE): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXE)
%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
clean:
	rm -f *.o *.mod                             
