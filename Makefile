DEBUGOPT = -Wall -g -Og -std=c11
FLINTOPT = -lm -lflint -lmpfr -lgmp -lpthread -I /usr/local/include/flint/
OUT = binding/main
TARGET = binding/main.c
CC = gcc
current_dir = $(shell pwd)
parent_dir = $(shell dirname "$(current_dir)")

compile:
	$(CC) $(TARGET) -o $(OUT) $(FLINTOPT) -O2 -std=c11

debug:
	$(CC) $(TARGET) -o $(OUT) $(DEBUGOPT) $(FLINTOPT)
clean:
	rm -f OUT

test-coset:
	sage -c 'sys.path.append("$(parent_dir)"); import e8theta_degree3.tests.hecke_module_test'
test-gen-gauss-sum:
	sage -c 'sys.path.append("$(parent_dir)"); import e8theta_degree3.tests.gen_gauss_sum_test'
test-minkowski:
	sage -c 'sys.path.append("$(parent_dir)"); import e8theta_degree3.tests.minkowski_reduction_test'

.PHONY: compile clean debug
