#ifndef COLOR_H
#define COLOR_H

class Color {
public:
	float r;
	float g;
	float b;

	Color ();
	Color (Color* col);
	Color (float r, float g, float b);

	Color operator+(Color c) const;
	Color operator*(Color c) const;
	Color operator*(float s) const;
};

#endif