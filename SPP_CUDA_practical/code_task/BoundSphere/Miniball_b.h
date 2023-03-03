#ifndef MINIBALL_B_H
#define MINIBALL_B_H

#include "PointN.h"
#include <cmath>

inline double mb_sqr (double r) {return r*r;}

class Miniball_b {
private:
    
  // data members
  int                 d;      // dimension
  int                 m, s;   // size and number of support points
  double*             q0;
    
  double*             z;
  double*             f;
  double**            v;
  double**            a;
    
  double**            c;
  double*             sqr_r;
    
  double*             current_c;      // refers to some c[j]
  double              current_sqr_r; 
 
  Miniball_b (const Miniball_b& mb);
  Miniball_b& operator=(const Miniball_b& mb); 

public:
  Miniball_b(int dim); 


  ~Miniball_b(); 
    
  // access
  const double*       center() const;
  double              squared_radius() const;
  int                 size() const;
  int                 support_size() const;
  double              excess (const PointN& p) const;
    
  // modification
  void                reset(); // generates empty sphere with m=s=0
  bool                push (const PointN& p);
  void                pop ();
    
  // checking
  double              slack() const;
};

#endif