#ifndef _MANDELBROT_H
#define _MANDELBROT_H

#include "image_distributed.h"

/*--- Type definitions -----------------------------------------------------*/

/**
 * This structure is used to pass the set of required parameters to the
 * mandelbrot() call.
 */
typedef struct {
  /* Input: section of the complex plane */
  double xmin; /**< Lower bound in complex plane (real part) */
  double ymin; /**< Lower bound in complex plane (imag. part) */
  double xmax; /**< Upper bound in complex plane (real part) */
  double ymax; /**< Upper bound in complex plane (imag. part) */
  int maxiter; /**< Maximum number of iterations */

  /* Input: image size & offsets */
  int columns; /**< Number of pixels to draw in x direction */
  int rows;    /**< Number of pixels to draw in y direction */

  // assigned to this process:
  int from; // inclusive
  int to;   // exclusive

  /* Output: image */
  image_t *image; /**< Pointer to image data structure */
} mandel_t;

/*--- Function prototypes --------------------------------------------------*/

void *mandelbrot(mandel_t *data);

#endif /* !_MANDELBROT_H */
