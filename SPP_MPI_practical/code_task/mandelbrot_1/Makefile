CC = mpicc
CFLAGS = -Wall -Wextra -O2 -fopenmp -g
LDLIBS = -lm

all : mandel


clean :
	rm -f mandel *.o

mandel: image_distributed.o main.o mandelbrot.o utility.o
	$(CC) $(CFLAGS) -o mandel image_distributed.o main.o mandelbrot.o utility.o $(LDLIBS)

image.o : image_distributed.c image_distributed.h
	$(CC) $(CFLAGS) -c image_distributed.c

main.o : main.c image_distributed.h mandelbrot.h
	$(CC) $(CFLAGS) -c main.c

mandelbrot.o : mandelbrot.c mandelbrot.h image_distributed.h utility.h
	$(CC) $(CFLAGS) -c mandelbrot.c

utility.o : utility.c utility.h image_distributed.h
	$(CC) $(CFLAGS) -c utility.c
