FC = ifort
#FCFLAGS = -g -fbounds-check 
FCFLAGS = -O3
#FCFLAGS += -I/
LIBS = -L/opt/lapack/lib -llapack -lblas -mkl

PROGRAMS = solver

solver.o: parameters.o 
solver: parameters.o 

%: %.o 
	$(FC) $(FCFLAGS) -o $@ $^ $(LIBS)

%.o : %.f90 
	$(FC) $(FCFLAGS) -c $< $(LIBS)

all:$(PROGRAMS)

.PHONY: clean

clean: 
	rm *.o *.mod
