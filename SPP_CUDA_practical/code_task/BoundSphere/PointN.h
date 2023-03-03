#ifndef POINTN_H
#define POINTN_H

#include <cassert>

class PointN {
private:
  int d; 
  double* coord;
   
public:
  // default
  PointN(int dim)
    : d (dim), coord(new double[dim])
  {}

  ~PointN ()
  {
    delete[] coord;
  }
   
  // copy from PointN
  PointN (const PointN& p)
    : d (p.dim()), coord(new double[p.dim()])
  {
    for (int i=0; i<d; ++i)
      coord[i] = p.coord[i];
  }
   
  // copy from double*
  PointN (int dim, const double* p)
    : d (dim), coord(new double[dim])
  {
    for (int i=0; i<d; ++i)
      coord[i] = p[i];
  }
   
  // assignment
  PointN& operator = (const PointN& p)
  {
    assert (d == p.dim());
    if (this != &p)
      for (int i=0; i<d; ++i)
	coord[i] = p.coord[i];
    return *this;
  }
   
  // coordinate access
  double& operator [] (int i)
  {
    return coord[i];
  }
  const double& operator [] (int i) const
  {
    return coord[i];
  }
  const double* begin() const
  {
    return coord;
  }
  const double* end() const
  {
    return coord+d;
  }

  // dimension access
  int dim() const
  {
    return d;
  }
};

#endif