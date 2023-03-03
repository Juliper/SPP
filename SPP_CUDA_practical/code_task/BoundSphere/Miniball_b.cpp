#include "Miniball_b.h"

Miniball_b::Miniball_b(int dim) : d(dim) {
  q0 = new double[d];
  z  = new double[d + 1];
  f  = new double[d + 1];
  v  = new double*[d + 1];
  for (int i = 0; i < d + 1; ++i)
    v[i] = new double[d];
  a = new double*[d + 1];
  for (int i = 0; i < d + 1; ++i)
    a[i] = new double[d];
  c = new double*[d + 1];
  for (int i = 0; i < d + 1; ++i)
    c[i] = new double[d];
  sqr_r = new double[d + 1];
  reset();
}

Miniball_b::~Miniball_b() {
  delete[] sqr_r;
  for (int i = 0; i < d + 1; ++i)
    delete[] c[i];
  delete[] c;
  for (int i = 0; i < d + 1; ++i)
    delete[] a[i];
  delete[] a;
  for (int i = 0; i < d + 1; ++i)
    delete[] v[i];
  delete[] v;
  delete[] f;
  delete[] z;
  delete[] q0;
}

const double* Miniball_b::center() const {
  return current_c;
}

double Miniball_b::squared_radius() const {
  return current_sqr_r;
}

int Miniball_b::size() const {
  return m;
}

int Miniball_b::support_size() const {
  return s;
}

double Miniball_b::excess(const PointN& p) const {
  double e = -current_sqr_r;
  for (int k = 0; k < d; ++k)
    e += mb_sqr(p[k] - current_c[k]);
  return e;
}

void Miniball_b::reset() {
  m = s = 0;
  // we misuse c[0] for the center of the empty sphere
  for (int j = 0; j < d; ++j)
    c[0][j] = 0;
  current_c     = c[0];
  current_sqr_r = -1;
}

void Miniball_b::pop() {
  --m;
}

bool Miniball_b::push(const PointN& p) {
  int i, j;
  double eps = 1e-32;
  if (m == 0) {
    for (i = 0; i < d; ++i)
      q0[i] = p[i];
    for (i = 0; i < d; ++i)
      c[0][i] = q0[i];
    sqr_r[0] = 0;
  } else {
    // set v_m to Q_m
    for (i = 0; i < d; ++i)
      v[m][i] = p[i] - q0[i];

    // compute the a_{m,i}, i< m
    for (i = 1; i < m; ++i) {
      a[m][i] = 0;
      for (j = 0; j < d; ++j)
        a[m][i] += v[i][j] * v[m][j];
      a[m][i] *= (2 / z[i]);
    }

    // update v_m to Q_m-\bar{Q}_m
    for (i = 1; i < m; ++i) {
      for (j = 0; j < d; ++j)
        v[m][j] -= a[m][i] * v[i][j];
    }

    // compute z_m
    z[m] = 0;
    for (j = 0; j < d; ++j)
      z[m] += mb_sqr(v[m][j]);
    z[m] *= 2;

    // reject push if z_m too small
    if (z[m] < eps * current_sqr_r) {
      return false;
    }

    // update c, sqr_r
    double e = -sqr_r[m - 1];
    for (i = 0; i < d; ++i)
      e += mb_sqr(p[i] - c[m - 1][i]);
    f[m] = e / z[m];

    for (i = 0; i < d; ++i)
      c[m][i] = c[m - 1][i] + f[m] * v[m][i];
    sqr_r[m] = sqr_r[m - 1] + e * f[m] / 2;
  }
  current_c     = c[m];
  current_sqr_r = sqr_r[m];
  s             = ++m;
  return true;
}

double Miniball_b::slack() const {
  double* l    = new double[d + 1];
  double min_l = 0;
  l[0]         = 1;
  for (int i = s - 1; i > 0; --i) {
    l[i] = f[i];
    for (int k = s - 1; k > i; --k)
      l[i] -= a[k][i] * l[k];
    if (l[i] < min_l)
      min_l = l[i];
    l[0] -= l[i];
  }
  if (l[0] < min_l)
    min_l = l[0];
  delete[] l;
  return ((min_l < 0) ? -min_l : 0);
}