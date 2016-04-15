DEBUGOPT = -Wall -g -Og -std=c11
FLINTOPT = -lm -lflint -lmpfr -lgmp -lpthread -I /usr/local/include/flint/
OUT = binding/main
TARGET = binding/main.c
CC = gcc
compile:
	$(CC) $(TARGET) -o $(OUT) $(FLINTOPT) -O2 -std=c11

debug:
	$(CC) $(TARGET) -o $(OUT) $(DEBUGOPT) $(FLINTOPT)
clean:
	rm -f OUT

.PHONY: compile clean debug
