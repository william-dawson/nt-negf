FC = mpif90
FCFLAGS = -fopenmp -fcheck=bounds
INCLUDES = -I./external/include -I./external/include/modfiles
LDFLAGS = -L./external/lib
LDLIBS = -lNTPoly -llapack -lblas -lnegf -lmpifx

MATS = example2/data-sw2/hamiltonian_sparse.mtx \
       example2/data-sw2/overlap_sparse.mtx \
       example2/data-sw2/density_kernel_sparse.mtx

EXEC = driv

.PHONY: all test1 test2 clean

all: $(EXEC)

$(EXEC): src/prototype.f90
	$(FC) $(FCFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS) $(LDLIBS)

test1: $(EXEC)
	mpirun -np 1 ./$< $(MATS) order.txt out.mtx

test2: $(EXEC)
	mpirun -np 2 ./$< $(MATS) order.txt out.mtx

clean:
	rm -f $(EXEC)
