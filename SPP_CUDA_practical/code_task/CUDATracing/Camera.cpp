#include "Camera.h"

#include "Constants.h"
#include "Log.h"

Camera::Camera(Point* pos, float cameraFOVwidth, float cameraFOVheight, int imageWidth, int imageHeight) {
  this->pos             = pos;
  this->cameraFOVwidth  = cameraFOVwidth;
  this->cameraFOVheight = cameraFOVheight;
  this->imageWidth      = imageWidth;
  this->imageHeight     = imageHeight;

  pixelWidth  = cameraFOVwidth / imageWidth;
  pixelHeight = -pixelWidth;
}

Camera::Camera() {
}

Ray Camera::getCamRay(int x, int y) {
  Direction dir = Direction(x * pixelWidth - 0.5f * imageWidth * pixelWidth,
                            y * pixelHeight - 0.5f * imageHeight * pixelHeight, -1.0f);
  return Ray(pos, &dir);
}

Ray3D Camera::getCamRay3D(int x, int y) {
  Direction dir = Direction(x * pixelWidth - 0.5f * imageWidth * pixelWidth,
                            y * pixelHeight - 0.5f * imageHeight * pixelHeight, -1.0);

  dir.normalize();
  return Ray3D(pos, dir);
}
