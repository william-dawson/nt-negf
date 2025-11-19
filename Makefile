FC = mpif90
FCFLAGS = -fopenmp -fcheck=bounds
INCLUDES = -I./external/include -I./external/include/modfiles
LDFLAGS = -L./external/lib
LDLIBS = -lNTPoly -llapack -lblas -lnegf -lmpifx

MATS = example/data-sw/hamiltonian_sparse.mtx \
       example/data-sw/overlap_sparse.mtx \
       example/data-sw/density_kernel_sparse.mtx \
       example/data-sw/order.txt

EXEC = driv

.PHONY: all test1 test2 clean

all: $(EXEC)

$(EXEC): src/prototype.f90
	$(FC) $(FCFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS) $(LDLIBS)

test1: $(EXEC)
	mpirun -np 1 ./$< $(MATS) out.mtx

test2: $(EXEC)
	mpirun -np 2 ./$< $(MATS) out.mtx

clean:
	rm -f $(EXEC)
