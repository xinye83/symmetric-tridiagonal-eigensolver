FC    = ifort
MPIFC = mpiifort

MKL= -Wl,--start-group $(MKLLIB)/libmkl_intel_ilp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a -Wl,--end-group

FLAG =-i8 -traceback -check bounds -qopenmp -lpthread -lm

MKLLIB= $(MKLROOT)/lib/intel64
MKLINC= $(MKLROOT)/include


all: main

main: mmio.o stack_mod.o subroutines.o main.o 
	$(MPIFC) -o $@ $^ -L$(MKLLIB) -I$(MKLINC) $(FLAG) $(MKL)

mmio.o: mmio.f
	$(MPIFC) -c $<

main.o: main.f90 mmio.o stack_mod.o
	$(MPIFC) -c $< -I$(MKLINC) $(FLAG)

stack_mod.o: stack_mod.f90
	$(MPIFC) -c $< -I$(MKLINC) $(FLAG)

subroutines.o: subroutines.f90
	$(MPIFC) -c $< -I$(MKLINC) $(FLAG)


clean:
	rm -f main *.o *.mod
