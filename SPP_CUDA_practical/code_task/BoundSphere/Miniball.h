#ifndef MINIBALL_H
#define MINIBALL_H

#include "Miniball_b.h"
#include "PointN.h"

#include <cassert>
#include <list>
#include <cmath>

class Miniball {
public:
  // types
  typedef std::list<PointN>::iterator         It;
  typedef std::list<PointN>::const_iterator   Cit;
    
private:
  // data members
  int              d;            // dimension
  std::list<PointN> L;            // internal point set
  Miniball_b       B;            // the current ball
  It               support_end;  // past-the-end iterator of support set
    
  // private methods
  void        mtf_mb (It k);
  void        pivot_mb (It k);
  void        move_to_front (It j);
  double      max_excess (It t, It i, It& pivot) const;  

public:
  // creates an empty ball
  Miniball(int dim) : d(dim), B(dim) {}

  // copies p to the internal point set
  void        check_in (const PointN& p);

  // builds the smallest enclosing ball of the internal point set
  void        build ();
    
  // returns center of the ball (undefined if ball is empty)
  PointN       center() const;

  // returns squared_radius of the ball (-1 if ball is empty)
  double      squared_radius () const;

  // returns size of internal point set
  int         nr_points () const;

  // returns begin- and past-the-end iterators for internal point set
  Cit         points_begin () const;
  Cit         points_end () const;

  // returns size of support point set; this set has the following properties:
  // - there are at most d+1 support points, 
  // - all support points are on the boundary of the computed ball, and
  // - the smallest enclosing ball of the support point set equals the 
  //   smallest enclosing ball of the internal point set
  int         nr_support_points () const;

  // returns begin- and past-the-end iterators for internal point set
  Cit         support_points_begin () const;
  Cit         support_points_end () const;
    
  // assesses the quality of the computed ball. The return value is the
  // maximum squared distance of any support point or point outside the 
  // ball to the boundary of the ball, divided by the squared radius of
  // the ball. If everything went fine, this will be less than e-15 and
  // says that the computed ball approximately contains all the internal
  // points and has all the support points on the boundary.
  // 
  // The slack parameter that is set by the method says something about
  // whether the computed ball is really the *smallest* enclosing ball 
  // of the support points; if everything went fine, this value will be 0; 
  // a positive value may indicate that the ball is not smallest possible,
  // with the deviation from optimality growing with the slack
  double      accuracy (double& slack) const;

  // returns true if the accuracy is below the given tolerance and the
  // slack is 0
  bool        is_valid (double tolerance = 1e-15) const;
};

#endif // !MINIBALL_H