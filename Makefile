current_dir = $(shell pwd)
parent_dir = $(shell dirname "$(current_dir)")
DEBUGOPT = -Wall -g -Og -std=c11
PATHOPT = -L$(current_dir)/binding -I/usr/local/include/flint/ -I$(current_dir)/binding
LIBOPTBASE = -lm -lflint -lmpfr -lgmp -lpthread
OPT = -O2 -std=c11
SHARED = -shared -fPIC
CC = gcc

compile-e8vector:
	$(CC) binding/e8vectors.c -o binding/libe8vectors.so $(PATHOPT) $(OPT) $(LIBOPTBASE) $(SHARED)

compile-miyawaki-theta:
	$(CC) binding/miyawaki_theta.c -o binding/libmiyawaki_theta.so $(PATHOPT) $(OPT) -le8vectors $(LIBOPTBASE) $(SHARED)
compile:
	$(CC) binding/miyawaki_theta.c -o binding/miyawaki $(PATHOPT) $(OPT) -le8vectors $(LIBOPTBASE)

compile-cython: compile-e8vector compile-miyawaki-theta
	sage -c 'sh.eval("cd binding; python setup.py build_ext -i")'

debug:
	$(CC) $(TARGET) -o $(OUT) $(DEBUGOPT) $(PATHOPT) $(LIBOPT)
clean:
	rm -f OUT

test-coset:
	sage -c 'sys.path.append("$(parent_dir)"); import e8theta_degree3.tests.hecke_module_test'
test-gen-gauss-sum:
	sage -c 'sys.path.append("$(parent_dir)"); import e8theta_degree3.tests.gen_gauss_sum_test'
test-minkowski:
	sage -c 'sys.path.append("$(parent_dir)"); import e8theta_degree3.tests.minkowski_reduction_test'

.PHONY: compile clean debug
