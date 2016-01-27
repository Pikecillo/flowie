/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __BLOBS_HPP__
#define __BLOBS_HPP__

#include <vector>

#include "vector_types.hpp"

/*template <class S_, int dim_>
class BlobbySurface {
public:
  BlobbySurface(S_ thresh=0.0) : threshold(thresh) {}

  void push_back(const Vector<S_, dim_> center, S_ weight, S_ radius) {
    centers.push_back(center);
    weights.push_back(weight);
    radii.push_back(radius);
  }

  S_ eval(const Vector<S_, dim_> &point) const {
    S_ f = 0.0;

    for(unsigned int i = 0; i < centers.size(); i++) {
      S_ r = (point - centers[i]).length();
      f -= (function(r) / radii[i] * weights[i] - threshold);
    }

    return f;
  }

private:
  S_ function(S_ r) const {
    S_ fac = (1 - r * r);
    return fac * fac * fac;
  }

private:
  S_ threshold;
  std::vector<Vector<S_, dim_> > centers;
  std::vector<S_> weights;
  std::vector<S_> radii;
};

typedef BlobbySurface<float, 2> BlobbySurface2f;
typedef BlobbySurface<float, 3> BlobbySurface3f;

typedef BlobbySurface<double, 2> BlobbySurface2d;
typedef BlobbySurface<double, 2> BlobbySurface2d;*/

#endif
