#ifndef _UTILITY_H
#define _UTILITY_H

//#include "image_distributed.h"

/*--- Function prototypes --------------------------------------------------*/

/**
 * Data type for a single pixel, capable of storing it's color in the RGB
 * color model.
 */
typedef struct {
  unsigned char red;   /**< Red color component */
  unsigned char green; /**< Green color component */
  unsigned char blue;  /**< Blue color component */
  unsigned char pad;   /**< "Padding" byte for alignment */
} color_t;

double get_wtime();
color_t HSVtoRGB(double hue, double sat, double val);

#endif /* !_UTILITY_H */
