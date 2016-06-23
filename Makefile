current_dir = $(shell pwd)
parent_dir = $(shell dirname "$(current_dir)")
DEBUGOPT = -Wall -g -Og -std=c11
PATHOPT = -L$(current_dir)/binding -I/usr/local/include/flint/
LIBOPT = -le8vectors -lm -lflint -lmpfr -lgmp -lpthread
OUT = binding/main
TARGET = binding/main.c
CC = gcc

compile:
	$(CC) $(TARGET) -o $(OUT) $(PATHOPT) -O2 -std=c11 $(LIBOPT)

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
