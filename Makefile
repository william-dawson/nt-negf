FC = mpif90
FCFLAGS = -fopenmp -fcheck=all -fbacktrace -g
INCLUDES = -I./external/include -I./external/include/modfiles
LDFLAGS = -L./external/lib
LDLIBS = -lNTPoly -llapack -lblas -lnegf -lmpifx

MATS = example/data-sw/hamiltonian_sparse.mtx \
       example/data-sw/overlap_sparse.mtx \
       example/data-sw/density_kernel_sparse.mtx \
       example/data-sw/order.txt

EXEC = driv
EXEC_SCF = driv_scf
SCF_MODULE = scf_module.o

.PHONY: all test1 test2 scf clean

all: $(EXEC) $(EXEC_SCF)

$(EXEC): src/prototype.f90
	$(FC) $(FCFLAGS) $(INCLUDES) $< -o $@ $(LDFLAGS) $(LDLIBS)

$(SCF_MODULE): src/scf_module.f90
	$(FC) $(FCFLAGS) $(INCLUDES) -c $< -o $@ $(LDFLAGS)

$(EXEC_SCF): $(SCF_MODULE) src/scf_driver.f90
	$(FC) $(FCFLAGS) $(INCLUDES) src/scf_driver.f90 $(SCF_MODULE) -o $@ $(LDFLAGS) $(LDLIBS)

test1: $(EXEC)
	mpirun -np 1 ./$< $(MATS) out.mtx

test2: $(EXEC)
	mpirun -np 2 ./$< $(MATS) out.mtx

scf: $(EXEC_SCF)
	mpirun -np 1 ./$(EXEC_SCF) --data-dir ./example/data-sw --max-iter 15 --tolerance 1e-3

scf2: $(EXEC_SCF)
	mpirun -np 1 ./$(EXEC_SCF) --data-dir ./example/data-sw2 --max-iter 15 --tolerance 1e-3

clean:
	rm -f $(EXEC) $(EXEC_SCF) $(SCF_MODULE) *.mod
