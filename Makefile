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

.PHONY: compile clean debug
