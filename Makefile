current_dir = $(shell pwd)
parent_dir = $(shell dirname "$(current_dir)")
DEBUGOPT = -Wall -g -Og -std=c11
PATHOPT = -L$(current_dir)/binding/lib -I$(current_dir)/binding/
OPT = -O3 -std=c11
SHARED = -shared -fPIC
CC = gcc

compile-theta_vectors:
	$(CC) binding/int_vector_sort.c -o binding/lib/libint_vector_sort.so $(PATHOPT) $(OPT) $(SHARED)
	$(CC) binding/e8vectors.c -o binding/lib/libe8vectors.so $(PATHOPT) -lint_vector_sort $(OPT) $(SHARED)
	$(CC) binding/rank16_vectors.c -o binding/lib/librank16_vectors.so $(PATHOPT) -lint_vector_sort $(OPT) $(SHARED)

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
