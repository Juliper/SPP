#include <stdio.h>
#include <stdlib.h>

#include <assert.h>
#include <mpi.h>

#include "image_distributed.h"

/*--- Implementation -------------------------------------------------------*/

/**
 * Allocates and initializes an image data structure of the given
 * dimensions.
 *
 * @param  width   Image width in pixels
 * @param  height  Image height in pixels
 *
 * @return Pointer to image data structure if successful, NULL otherwise
 */
image_t *imageCreate(int global_width, int global_height, int local_width,
                     int local_height, int x_offset, int y_offset) {
  int x;
  image_t *image;

  /* Allocate image data structure */
  image = (image_t *)malloc(sizeof(image_t));
  if (!image) {
    fprintf(stderr, "Memory allocation error!\n");
    return NULL;
  }

  /* Allocate imaga data array */
  image->data = (color_t **)malloc(local_width * sizeof(color_t *));
  if (!image->data) {
    fprintf(stderr, "Memory allocation error!\n");
    free(image);
    return NULL;
  }
  for (x = 0; x < local_width; ++x) {
    image->data[x] = (color_t *)malloc(local_height * sizeof(color_t));
    if (!image->data[x]) {
      int i;

      fprintf(stderr, "Memory allocation error!\n");
      for (i = x - 1; i >= 0; --i)
        free(image->data[i]);
      free(image->data);
      free(image);
      return NULL;
    }
  }

  /* Set attributes */
  image->global_width = global_width;
  image->global_height = global_height;

  image->local_width = local_width;
  image->local_height = local_height;

  image->x_offset = x_offset;
  image->y_offset = y_offset;

  return image;
}

/**
 * Releases all resources occupied by the given image data structure.
 *
 * @param  image  Image data structure to be freed
 */
void imageFree(image_t *image) {
  int x;

  /* Free up resources */
  for (x = 0; x < image->local_width; ++x)
    free(image->data[x]);
  free(image->data);
  free(image);
}

/**
 * Set the color of a single pixel.
 *
 * @param  image  Image data structure
 * @param  x      X coordinate
 * @param  y      Y coordinate
 * @param  color  Color value
 */
void imageSetPixel(image_t *image, int x, int y, color_t color) {
  // assert that acess is to a local px
  assert(x - image->x_offset >= 0 && x - image->x_offset < image->local_width);
  assert(y - image->y_offset >= 0 && y - image->y_offset < image->local_height);
  image->data[x - image->x_offset][y - image->y_offset] = color;
}

/**
 * Writes the given image to a PPM file with the provided name.
 *
 * @param  image     Image data structure
 * @param  filename  Name of output file
 */
void imageSave(image_t *image, const char *filename) {
  int x;
  int y;

  int rank, numprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  MPI_File file;

  if (rank == 0) { // only rank 0 writes the header
    FILE *fp;
    /* Open output file */
    fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Could not create output file \"%s\"!\n", filename);
      return;
    }

    /* Write PPM header */
    fprintf(fp, "P6\n%d %d\n255\n", image->global_width, image->global_height);
    fclose(fp);
  }

  // barrier to ensure that file was properly closed on rank 0 before any
  // process may open it
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_open(MPI_COMM_WORLD, filename,
                MPI_MODE_WRONLY | MPI_MODE_EXCL | MPI_MODE_APPEND,
                MPI_INFO_NULL, &file);

  // offset where the local part begins (assuming a row wise distribution of the
  // img)
  int file_offset = image->y_offset * image->global_width * 3;

  // goto the part assigned to this process
  MPI_File_seek(file, file_offset, MPI_SEEK_END);
  // seek after end: after header

  /* Write PPM data */
  for (y = 0; y < image->local_height; ++y)
    for (x = 0; x < image->local_width; ++x) {

      MPI_File_write(file, &(image->data[x][y]), 3, MPI_CHAR,
                     MPI_STATUS_IGNORE);
    }

  /* Close output file */
  MPI_File_close(&file);
}
