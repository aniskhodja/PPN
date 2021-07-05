CC = mpicc
CFLAGS = -Wall -lm -fopenmp -g -O3 -lm
all:./run
./run:main.c
	$(CC) $< -o $@ $(CFLAGS)
clean:
	rm -f *.o
	rm  run