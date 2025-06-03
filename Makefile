FC = mpif90
FCFLAGS = -fopenmp -fcheck=bounds
INCLUDES = -I./include -I./include/modfiles
LDFLAGS = -L./lib
LDLIBS = -lNTPoly -llapack -lblas -lnegf -lmpifx

MATS = example/data-sw/hamiltonian_sparse.mtx \
       example/data-sw/overlap_sparse.mtx \
       example/data-sw/density_kernel_sparse.mtx

EXEC = driv

.PHONY: all test1 test2 clean

all: $(EXEC)

$(EXEC): prototype.f90
	$(FC) $(FCFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS) $(LDLIBS)

test1: $(EXEC)
	mpirun -np 1 ./$< $(MATS) order.txt out.mtx

test2: $(EXEC)
	mpirun -np 2 ./$< $(MATS) order.txt out.mtx

clean:
	rm -f $(EXEC)
