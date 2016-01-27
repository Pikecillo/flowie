/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __INTERPOLANT_HPP__
#define __INTERPOLANT_HPP__

#include "vector_types.hpp"

/*
 * Linearly interpolates the value of a 1D function f(x) = y,
 * given the points (x0, y0) and (x1, y1)
 */
template <class D_, class R_>
inline R_ linear_interp(const D_ &x, const D_ &x0, const D_ &x1,
                        const R_ &y0, const R_ &y1) {
  D_ t = (x - x0) / (x1 - x0);
  return y0 + (y1 - y0) * t;
}

/*
 * Linearly interpolate the value of a 2D function f(x1, x2) = y,
 * given four points on a rectangular grid:
 * (x1[0], x2[0], y[0]), (x1[1], x2[0], y[1]),
 * (x1[0], x2[1], y[2]), (x1[1], x2[1], y[3])
 */
template <class D_, class R_>
inline R_ bilinear_interp(const Vector<D_, 2> &x,
                          const Vector<D_, 2> &x1,
                          const Vector<D_, 2> &x2,
                          R_ y[]) {
  R_ y01 = linear_interp(x[0], x1[0], x1[1], y[0], y[1]);
  R_ y23 = linear_interp(x[0], x1[0], x1[1], y[2], y[3]);
  return linear_interp(x[1], x2[0], x2[1], y01, y23);
}

/*
 * Linearly interpolate the value of a 3D function f(x1, x2, x3) = y,
 * given eight points on a rectangular grid:
 * (x1[0], x2[0], x3[0], y[0]), (x1[1], x2[0], x3[0], y[1]),
 * (x1[0], x2[1], x3[0], y[2]), (x1[1], x2[1], x3[0], y[3]),
 * (x1[0], x2[0], x3[1], y[4]), (x1[1], x2[0], x3[1], y[5]),
 * (x1[0], x2[1], x3[1], y[6]), (x1[1], x2[1], x3[1], y[7]),
 */
template <class D_, class R_>
R_ trilinear_interp(const Vector<D_, 3> &x,
                    const Vector<D_, 2> &x1,
                    const Vector<D_, 2> &x2,
                    const Vector<D_, 2> &x3,
                    R_ y[]) {
  R_ y01 = linear_interp(x[0], x1[0], x1[1], y[0], y[1]);
  R_ y23 = linear_interp(x[0], x1[0], x1[1], y[2], y[3]);
  R_ y45 = linear_interp(x[0], x1[0], x1[1], y[4], y[5]);
  R_ y67 = linear_interp(x[0], x1[0], x1[1], y[6], y[7]);
  R_ y_[] = { y01, y23, y45, y67 };
  Vector<D_, 2> x_(x[1], x[2]);
  return bilinear_interp(x_, x2, x3, y_);
}

#endif
