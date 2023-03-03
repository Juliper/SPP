#include "Miniball.h"

void Miniball::check_in(const PointN& p) {
  assert(d == p.dim());
  L.push_back(p);
}

void Miniball::build() {
  B.reset();
  support_end = L.begin();
  pivot_mb(L.end());
}

void Miniball::mtf_mb(It i) {
  support_end = L.begin();
  if ((B.size()) == d + 1)
    return;
  for (It k = L.begin(); k != i;) {
    It j = k++;
    if (B.excess(*j) > 0) {
      if (B.push(*j)) {
        mtf_mb(j);
        B.pop();
        move_to_front(j);
      }
    }
  }
}

void Miniball::move_to_front(It j) {
  if (support_end == j)
    support_end++;
  L.splice(L.begin(), L, j);
}

void Miniball::pivot_mb(It i) {
  It t = ++L.begin();
  mtf_mb(t);
  double max_e, old_sqr_r = -1;
  do {
    It pivot;
    max_e = max_excess(t, i, pivot);
    if (max_e > 0) {
      t = support_end;
      if (t == pivot)
        ++t;
      old_sqr_r = B.squared_radius();
      B.push(*pivot);
      mtf_mb(support_end);
      B.pop();
      move_to_front(pivot);
    }
  } while ((max_e > 0) && (B.squared_radius() > old_sqr_r));
}

double Miniball::max_excess(It t, It i, It& pivot) const {
  const double *c = B.center(), sqr_r = B.squared_radius();
  double e, max_e                     = 0;
  for (It k = t; k != i; ++k) {
    const double* p = (*k).begin();
    e               = -sqr_r;
    for (int j = 0; j < d; ++j)
      e += mb_sqr(p[j] - c[j]);
    if (e > max_e) {
      max_e = e;
      pivot = k;
    }
  }
  return max_e;
}

PointN Miniball::center() const {
  return PointN(d, B.center());
}

double Miniball::squared_radius() const {
  return B.squared_radius();
}

int Miniball::nr_points() const {
  return (int)L.size();
}

Miniball::Cit Miniball::points_begin() const {
  return L.begin();
}

Miniball::Cit Miniball::points_end() const {
  return L.end();
}

int Miniball::nr_support_points() const {
  return B.support_size();
}

Miniball::Cit Miniball::support_points_begin() const {
  return L.begin();
}

Miniball::Cit Miniball::support_points_end() const {
  return support_end;
}

double Miniball::accuracy(double& slack) const {
  double e, max_e = 0;
  int n_supp = 0;
  Cit i;
  for (i = L.begin(); i != support_end; ++i, ++n_supp)
    if ((e = std::abs(B.excess(*i))) > max_e)
      max_e = e;

  // you've found a non-numerical problem if the following ever fails
  assert(n_supp == nr_support_points());

  for (i = support_end; i != L.end(); ++i)
    if ((e = B.excess(*i)) > max_e)
      max_e = e;

  slack = B.slack();
  return (max_e / squared_radius());
}

bool Miniball::is_valid(double tolerance) const {
  double slack;
  return ((accuracy(slack) < tolerance) && (slack == 0));
}