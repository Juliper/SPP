#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

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
void master_main(mandel_t *data);

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

  char *filename = "output.ppm";

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

  if (rank == 0) { // only rank 0 writes the header
    FILE *fp;
    /* Open output file */
    fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Could not create output file \"%s\"!\n", filename);
      return EXIT_FAILURE;
    }

    /* Write PPM header */
    fprintf(fp, "P6\n%d %d\n255\n", IMG_WIDTH, IMG_HEIGHT);
    fclose(fp);
  }

  // barrier to ensure that file was properly closed on rank 0 before any
  // process may open it
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_open(MPI_COMM_WORLD, filename,
                MPI_MODE_WRONLY | MPI_MODE_EXCL | MPI_MODE_APPEND,
                MPI_INFO_NULL, &(data->file));

  if (rank == 0) { // master
    master_main(data);
  } else { // worker
    mandelbrot(data);
  }

  MPI_File_close(&(data->file));
  free(data);

  MPI_Finalize();
  return EXIT_SUCCESS;
}

void master_main(mandel_t *data) {

  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  int buffer = 0;

  MPI_Status status;

  for (int row = 0; row < data->rows; row++) {
    MPI_Recv(&buffer, 1, MPI_INT, MPI_ANY_SOURCE, MESSAGE_TAG, MPI_COMM_WORLD,
             &status);
    int from = status.MPI_SOURCE;
    // tell him which row to calculate next
    buffer = row;
    MPI_Send(&buffer, 1, MPI_INT, from, MESSAGE_TAG, MPI_COMM_WORLD);
  }

  // numprocs -1 iterations as one do not message itself
  for (int i = 0; i < numprocs - 1; i++) {
    MPI_Recv(&buffer, 1, MPI_INT, MPI_ANY_SOURCE, MESSAGE_TAG, MPI_COMM_WORLD,
             &status);
    int from = status.MPI_SOURCE;
    // tell him to stop
    buffer = -1;
    MPI_Send(&buffer, 1, MPI_INT, from, MESSAGE_TAG, MPI_COMM_WORLD);
  }
}
