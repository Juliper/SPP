#ifndef _IMAGE_DISTRIBUTED_H
#define _IMAGE_DISTRIBUTED_H

/*--- Type definitions -----------------------------------------------------*/

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

/**
 * Data type for an RGB image of a specific size.
 */
typedef struct {
  int global_width;  /**< Image width in pixel */
  int global_height; /**< Image height in pixel */
  int local_width;
  int local_height;
  int x_offset;
  int y_offset;
  color_t **data; /**< Image data (2D array of pixel values) */
} image_t;

/*--- Function prototypes --------------------------------------------------*/

image_t *imageCreate(int global_width, int global_height, int local_width,
                     int local_height, int x_offset, int y_offset);
void imageFree(image_t *image);
void imageSetPixel(image_t *image, int x, int y, color_t color);
void imageSave(image_t *image, const char *filename);

#endif /* !_IMAGE_DISTRIBUTED_H */
