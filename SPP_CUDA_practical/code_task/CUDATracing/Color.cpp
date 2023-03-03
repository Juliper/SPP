#include "Color.h"

#include <math.h>

Color::Color() {
  this->r = 0;
  this->g = 0;
  this->b = 0;
}

Color::Color(Color* col) {
  this->r = col->r;
  this->g = col->g;
  this->b = col->b;
}

Color::Color(float r, float g, float b) {
  if (r > 1.0f)
    this->r = 1.0f;
  else if (r < 0.0f)
    this->r = 0.0f;
  else
    this->r = r;

  if (g > 1.0f)
    this->g = 1.0f;
  else if (g < 0.0f)
    this->g = 0.0f;
  else
    this->g = g;

  if (b > 1.0f)
    this->b = 1.0f;
  else if (b < 0.0f)
    this->b = 0.0f;
  else
    this->b = b;
}

Color Color::operator*(Color col) const {
  float r = this->r * col.r;
  if (r > 1.0f)
    r = 1.0f;

  float b = this->b * col.b;
  if (b > 1.0f)
    b = 1.0f;

  float g = this->g * col.g;
  if (g > 1.0f)
    g = 1.0f;

  return Color(r, g, b);
}

Color Color::operator*(float s) const {
  float r = this->r * s;
  if (r > 1.0f)
    r = 1.0f;

  float b = this->b * s;
  if (b > 1.0f)
    b = 1.0f;

  float g = this->g * s;
  if (g > 1.0f)
    g = 1.0f;

  return Color(r, g, b);
}

Color Color::operator+(Color col) const {
  float r = this->r + col.r;
  if (r > 1.0f)
    r = 1.0f;

  float b = this->b + col.b;
  if (b > 1.0f)
    b = 1.0f;

  float g = this->g + col.g;
  if (g > 1.0f)
    g = 1.0f;

  return Color(r, g, b);
}
