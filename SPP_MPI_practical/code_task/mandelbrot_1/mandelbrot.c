#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mandelbrot.h"
#include "utility.h"

/*--- Implementation -------------------------------------------------------*/

/**
 * Calculates an image of the mandelbrot set for the parameters given in
 * @p data (see description of mandel_t for details). This function takes
 * ownership of the @p data provided and releases the data structure after
 * finishing the calculation. Also, this function prints the wall-clock
 * time required to do the calculations.
 *
 * @param  data  Mandelbrot parameters
 *
 * @return Always NULL
 */
void *mandelbrot(mandel_t *data) {
  int x;
  int y;
  double dx;
  double dy;
  double start_time;
  double end_time;

  /* Time measurement */
  start_time = get_wtime();

  /* Initialization */
  dx = (data->xmax - data->xmin) / data->columns;
  dy = (data->ymax - data->ymin) / data->rows;

  /* Iterate over all rows */
  // meaning iterate over space for this process only
  for (y = data->from; y < data->to; ++y) {
    double c_imag = data->ymin + (y * dy);

    /* Iterate over all columns */
    for (x = 0; x < data->columns; ++x) {
      color_t color;
      double c_real = data->xmin + (x * dx);
      int iter = 0;
      double z_real = 0.0;
      double z_imag = 0.0;
      double z_norm = 0.0;

      /* Check whether recursive Mandelbrot equation remains bounded */
      // The actual calculation
      while (z_norm < 4.0 && iter < data->maxiter) {
        double z2_real = (z_real * z_real) - (z_imag * z_imag);
        double z2_imag = (z_real * z_imag) + (z_imag * z_real);

        z_real = z2_real + c_real;
        z_imag = z2_imag + c_imag;
        z_norm = (z_real * z_real) + (z_imag * z_imag);

        ++iter;
      }

      /* Bounded => black */
      if (iter == data->maxiter) {
        color.red = 0;
        color.green = 0;
        color.blue = 0;
        color.pad = 0;
      }
      /* Unbounded => compute nice color */
      else {
        color = HSVtoRGB(sqrt((double)iter / data->maxiter), 0.8, 0.8);
      }

      imageSetPixel(data->image, x, y, color);
    }
  }

  /* Time measurement */
  end_time = get_wtime();
  printf("Calculation time: %2.6f seconds\n", end_time - start_time);

  return NULL;
}
