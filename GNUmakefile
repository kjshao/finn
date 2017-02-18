FC = ifort
FFLAGS = -O3 -openmp -warn -module mod
LIB = -L/work1/soft/intel2015/composer_xe_2015.2.164/mkl/lib/intel64 -lmkl_rt
SRC = $(wildcard src/*)
OBJf90 = $(patsubst src/%.f90, %.o, $(filter %.f90, $(SRC)))
OBJmod = $(patsubst src/%.F90, %.m.o, $(filter %.F90, $(SRC)))
OBJcud = $(patsubst src/%.cu, %.cu.o, $(filter %.cu, $(SRC)))
EXE = bin/invariants

$(EXE): $(patsubst %.o, obj/%.o, $(OBJmod) $(OBJf90))
	$(FC) $(FFLAGS) $(LIB) \
	$+ -o $@
obj/%.o : src/%.f90
	$(FC) -c $< $(FFLAGS) $(LIB) -o $@
obj/%.m.o : src/%.F90
	$(FC) -c $< $(FFLAGS) $(LIB) -o $@
clean:
	rm -f obj/*.o mod/*.mod $(EXE)
