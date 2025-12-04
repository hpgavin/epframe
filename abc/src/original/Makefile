# Makefile for EPFRAME - Elastic-Plastic Frame Analysis
# Usage: make          (compile)
#        make test     (run example)
#        make clean    (remove executables)

CC = gcc
CFLAGS = -O2 -Wall -lm
PYTHON = python3

# C version
epframe: epframe.c NRutil.c NRutil.h
	$(CC) -o epframe epframe.c NRutil.c $(CFLAGS)

# Test with example input file
test: epframe
	./epframe example_frame.dat output_c.txt
	@echo ""
	@echo "=== C Version Output ==="
	@head -30 output_c.txt
	@echo ""
	@echo "Python version:"
	$(PYTHON) epframe.py example_frame.dat output_py.txt
	@echo ""
	@echo "=== Python Version Output ==="
	@head -30 output_py.txt

# Create example input file
example_frame.dat:
	@echo "5" > example_frame.dat
	@echo "8   7   29000" >> example_frame.dat
	@echo "1     0     0   0   0   0" >> example_frame.dat
	@echo "2     0   168   1   1   1" >> example_frame.dat
	@echo "3   120   252   1   1   1" >> example_frame.dat
	@echo "4   216   252   1   1   1" >> example_frame.dat
	@echo "5   312   252   1   1   1" >> example_frame.dat
	@echo "6   408   252   1   1   1" >> example_frame.dat
	@echo "7   528   168   1   1   1" >> example_frame.dat
	@echo "8   528     0   0   0   0" >> example_frame.dat
	@echo "1   1   2   954   19.7   4680" >> example_frame.dat
	@echo "2   2   3   954   19.7   4680" >> example_frame.dat
	@echo "3   3   4   954   19.7   4680" >> example_frame.dat
	@echo "4   4   5   954   19.7   4680" >> example_frame.dat
	@echo "5   5   6   954   19.7   4680" >> example_frame.dat
	@echo "6   6   7   954   19.7   4680" >> example_frame.dat
	@echo "7   7   8   954   19.7   4680" >> example_frame.dat
	@echo "5" >> example_frame.dat
	@echo "2   0.5   0   0" >> example_frame.dat
	@echo "3   0.25   -1   0" >> example_frame.dat
	@echo "4   0   -1   0" >> example_frame.dat
	@echo "5   0   -1   0" >> example_frame.dat
	@echo "6   0   -1   0" >> example_frame.dat
	@echo "Example input file created: example_frame.dat"

clean:
	rm -f epframe *.o output_c.txt output_py.txt

.PHONY: test clean example_frame.dat
