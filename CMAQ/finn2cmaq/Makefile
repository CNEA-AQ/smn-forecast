.SUFFIXES: .o .f90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm 
INC  = -I/usr/include
FFLAGS = -O2 -ffree-line-length-none
OBJS = finn2cmaq.o

.f90.o:
		${FC} ${FFLAGS} -c ${INC} $<

finn2cmaq: ${OBJS}
		${FC} -o $@ ${OBJS} ${LIBS} 

clean:
		rm -f *.o *.mod
