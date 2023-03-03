CC = mpicc
CFLAGS = -Wall -Wextra -O2 -fopenmp -g
LDLIBS = -lm

all : mandel


clean :
	rm -f mandel *.o

mandel: main.o mandelbrot.o utility.o
	$(CC) $(CFLAGS) -o mandel main.o mandelbrot.o utility.o $(LDLIBS)

main.o : main.c mandelbrot.h
	$(CC) $(CFLAGS) -c main.c

mandelbrot.o : mandelbrot.c mandelbrot.h utility.h
	$(CC) $(CFLAGS) -c mandelbrot.c

utility.o : utility.c utility.h
	$(CC) $(CFLAGS) -c utility.c
