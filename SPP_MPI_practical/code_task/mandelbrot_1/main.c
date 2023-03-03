#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "image_distributed.h"
#include "mandelbrot.h"

/** Width of output image in pixels */
#define IMG_WIDTH 4096

/** Height of output image in pixels */
#define IMG_HEIGHT 4096

/** Maximum number of iterations to perform */
#define MAX_ITER 5000

/**
 * Main program.
 */
int main() {
  MPI_Init(NULL, NULL);

  int rank, numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  /* Parameters */

  // teil der komplexen ebene, der betrachtet werden soll:  (globale werte)
  double xmin = -1.172657;
  double ymin = -0.296335;
  double xmax = -1.172643;
  double ymax = -0.296321;

  // use a row wise distribution of the img among processes

  int own_height = IMG_HEIGHT / numprocs;
  int offset = own_height * rank;
  if (rank == numprocs - 1) {
    own_height = IMG_HEIGHT - offset;
    // assign him all the remaining lines
  }
  // printf for debug:
  // printf("rank %d: %d lines starting with %d\n", rank, own_height, offset);

  /* Create image data structure */
  image_t *image =
      imageCreate(IMG_WIDTH, IMG_HEIGHT, IMG_WIDTH, own_height, 0, offset);
  if (!image) {
    fprintf(stderr, "Memory allocation error!\n");
    return EXIT_FAILURE;
  }

  /* Allocate mandelbrot data structure */
  mandel_t *data = (mandel_t *)malloc(sizeof(mandel_t));
  if (!data) {
    fprintf(stderr, "Memory allocation error!\n");
    return EXIT_FAILURE;
  }

  /* Initialize data */
  data->xmin = xmin;
  data->ymin = ymin;
  data->xmax = xmax;
  data->ymax = ymax;
  data->maxiter = MAX_ITER;
  data->columns = IMG_WIDTH;
  data->rows = IMG_HEIGHT;
  data->from = offset;
  data->to = offset + own_height;
  data->image = image;

  mandelbrot(data);

  free(data);

  /* Save the output image & free resources */
  imageSave(image, "output.ppm");
  imageFree(image);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
