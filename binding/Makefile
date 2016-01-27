DEBUGOPT = -Wall -g -O0
FLINTOPT = -lm -lflint -lmpfr -lgmp -lpthread -I /usr/local/include/flint/
OUT = main
TARGET = main.c
CC = gcc
compile:
	$(CC) $(TARGET) -o $(OUT) $(FLINTOPT) -Ofast

debug:
	$(CC) $(TARGET) -o $(OUT) $(DEBUGOPT) $(FLINTOPT)
clean:
	rm -f OUT

.PHONY: compile clean debug
