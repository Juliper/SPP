#include <stdio.h>
#include <sys/time.h>

#include "utility.h"

/*--- Implementation -------------------------------------------------------*/

/**
 * Returns a wall-clock time value in seconds since the Epoch (i.e.,
 * 1970-01-01 00:00:00 +0000).
 */
double get_wtime() {
  struct timeval tv;

  if (0 != gettimeofday(&tv, NULL))
    fprintf(stderr, "Error calling gettimeofday()\n");

  return tv.tv_sec + tv.tv_usec * 1e-6;
}

/**
 * Converts a color value given in the HSV color model into a color value in
 * the RGB color model. All values are expected to be in the interval [0,1].
 *
 * @param  hue  Hue value
 * @param  sat  Saturation value
 * @param  val  Value
 *
 * @return Corresponding color in RGB color model
 */
color_t HSVtoRGB(double hue, double sat, double val) {
  color_t result;
  unsigned int interval;
  double f;
  double p;
  double q;
  double t;
  double r;
  double g;
  double b;

  /* Clamp hue, saturation and value to interval [0,1] */
  if (hue < 0.0)
    hue = 0.0;
  if (hue > 1.0)
    hue = 1.0;
  if (sat < 0.0)
    sat = 0.0;
  if (sat > 1.0)
    sat = 1.0;
  if (val < 0.0)
    val = 0.0;
  if (val > 1.0)
    val = 1.0;

  /* Special case: hue = 0 */
  if (hue == 0.0) {
    result.red = 0.0;
    result.green = 0.0;
    result.blue = 0.0;

    return result;
  }

  /* Determine color interval */
  hue = hue * 6.0;
  if (hue == 6.0)
    hue = 0.0;
  interval = hue;

  /* Calculate intermediate values */
  f = hue - interval;
  p = val * (1 - sat);
  q = val * (1 - sat * f);
  t = val * (1 - sat * (1 - f));

  /* Determine RGB color value */
  switch (interval) {
  case 0:
    r = val;
    g = t;
    b = p;
    break;
  case 1:
    r = q;
    g = val;
    b = p;
    break;
  case 2:
    r = p;
    g = val;
    b = t;
    break;
  case 3:
    r = p;
    g = q;
    b = val;
    break;
  case 4:
    r = t;
    g = p;
    b = val;
    break;
  case 5:
    r = val;
    g = p;
    b = q;
    break;
  default:
    /* This should never happen... */
    r = 0.0;
    g = 0.0;
    b = 0.0;
    break;
  }

  /* Stretch color values from [0,1] to [0,255] */
  result.red = 255.0 * r;
  result.green = 255.0 * g;
  result.blue = 255.0 * b;
  result.pad = 0;

  return result;
}
