current_dir = $(shell pwd)
parent_dir = $(shell dirname "$(current_dir)")
DEBUGOPT = -Wall -g -Og -std=c11
PATHOPT = -L$(parent_dir)/lib -I$(current_dir) -I$(parent_dir)
LIBOPTBASE = -lmpir
OPT = -O3 -std=c11 -Wall -Wextra
SHARED = -shared -fPIC
CC = gcc

compile-c-lib:
	$(CC) theta14_13_5_c.c  -o $(parent_dir)/lib/libtheta14_13_5_c.so $(PATHOPT) $(OPT) \
	$(LIBOPTBASE) $(SHARED)

compile-cython: compile-c-lib
	sage -c 'sh.eval("python setup.py build_ext -i")'
