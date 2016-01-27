/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __RBF_HPP__
#define __RBF_HPP__

#include <assert.h>
#include <math.h>

#include <iostream>
#include <vector>

#include "vector_types.hpp"
#include "numerical/linear_solver.hpp"
#include "utils/reference.hpp"

template <class S_>
class RBFKernel {
public:
  virtual S_ eval(S_ r) = 0;
};

template <class S_>
class RBFBiharmonic : public RBFKernel<S_> {
public:
  virtual S_ eval(S_ r) {
    return r;
  }
};

template <class S_>
class RBFStam : public RBFKernel<S_> {
public:
  virtual S_ eval(S_ r) {
    S_ f = (1 - r * r);
    return f * f * f;
  }
};

template <class S_>
class RBFMultiquadric : public RBFKernel<S_> {
public:
  RBFMultiquadric(S_ s) : scale(s) {}

  virtual S_ eval(S_ r) {
    return sqrt(r * r + scale * scale);
  }

private:
  S_ scale;
};

template <class S_>
class RBFThinPlateSpline : public RBFKernel<S_> {
public:
  RBFThinPlateSpline(S_ s) : scale(s) {}

  virtual S_ eval(S_ r) {
    if(fabs(r) < 1E-10)
      return 0.0f;
    return r * r * log(r / scale);
  }

private:
  S_ scale;
};

template <class S_>
class RBFGaussian : public RBFKernel<S_> {
public:
  RBFGaussian(S_ s) : scale(s) {}

  virtual S_ eval(S_ r) {
    return exp(-0.5 * r * r / (scale * scale));
  }

private:
  S_ scale;
};

template <class S_, int dim_>
class RBF {
private:
  std::vector<Vector<S_, dim_> > xcenter; // Centers
  std::vector<S_> weights;
  std::vector<S_> ycenter;                // Values at centers
  RBFKernel<S_> *kernel;

public:
  RBF() : xcenter(), weights(), ycenter(), kernel(0) {}

  RBF(const RBF<S_, dim_> &other) {
    (*this) = other;
  }

  RBF(const std::vector<Vector<S_, dim_> > &x,
      const std::vector<S_> &y, RBFKernel<S_> *k)
    : xcenter(x), weights(), ycenter(y), kernel(k) {
    // Make sure there is something in the vectors
    assert(x.size() == y.size() && x.size() != 0);

    solveWeights();
  }

  void operator=(const RBF<S_, dim_> &other) {
    xcenter = other.xcenter;
    weights = other.weights;
    ycenter = other.ycenter;
    kernel = other.kernel;
  }

  S_ eval(const Vector<S_, dim_> &x) const {
    assert(kernel);

    S_ f = (S_)0;
    
    // Make sure the RBF has been initialized
    assert(weights.size());

    for(unsigned int i = 0; i < xcenter.size(); i++)
      f += weights[i] * kernel->eval((x - xcenter[i]).length());
  
    return f;
  }

  std::vector<Vector<S_, dim_> > getCenters() const { return xcenter; }
  std::vector<S_> getWeights() const { return weights; }

private:
  void solveWeights() {
    int n = ycenter.size();
    Matrixd y(n, 1), w(n, 1), A(n, n);

    for(int i = 0; i < n; i++) {
      y.at(i, 0) = ycenter[i];
      
      for(int j = 0; j < n; j++)
        A.at(i, j) = kernel->eval((xcenter[i] - xcenter[j]).length());
    }
    
    Solver::dgesv(A, y, w);
    
    for(int i = 0; i < n; i++)
      weights.push_back(w.at(i, 0));
  }
};

typedef RBF<float, 3> RBF3f;
typedef RBF<double, 3> RBF3d;

#endif
