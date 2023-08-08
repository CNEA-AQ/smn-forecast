.SUFFIXES: .o .f90

FC   = gfortran
LIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf -lnetcdff -lm
INC  = -I/usr/include
FFLAGS = -O2 -ffree-line-length-none -Wunused -x f95-cpp-input -DWGET #-g
OBJS = remap_test.o
EXE  = test.exe

.f90.o:
	${FC} ${FFLAGS} -c ${INC} $<

finn2cmaq: ${OBJS}
	${FC} -o ${EXE} ${OBJS} ${LIBS} 

clean:
	rm -f *.o *.mod

