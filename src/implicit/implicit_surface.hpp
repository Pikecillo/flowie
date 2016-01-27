/*
 * University of Houston
 * Mario Rincon-Nigro. May 2013.
 */

#ifndef __IMPLICIT_SURFACE_HPP__
#define __IMPLICIT_SURFACE_HPP__

#include <iostream>
#include <vector>

#include "vector_types.hpp"
#include "matrix_types.hpp"
#include "field_types.hpp"

template <class S_, int dim_>
class ImplicitSurface {
public:
  virtual S_ evaluate(const Vector<S_, dim_> &x) const = 0;
  virtual Vector<S_, dim_> gradient(const Vector<S_, dim_> &x) const = 0;
  virtual Mat<S_, dim_, dim_> hessian(const Vector<S_, dim_> &x) const = 0;
};

template <class S_, int dim_>
class DynImplicitSurface {
public:
  virtual S_ evaluate(const Vector<S_, dim_> &x, S_ t) const = 0;
  virtual Vector<S_, dim_> gradient(const Vector<S_, dim_> &x, S_ t) const = 0;
  virtual Mat<S_, dim_, dim_> hessian(const Vector<S_, dim_> &x, S_ t) const = 0;
  virtual S_ timeParDer(const Vector<S_, dim_> &x, S_ t) const = 0;
  virtual Vector<S_, dim_> timeParDerOfGradient(const Vector<S_, dim_> &x,
                                                S_ t) const = 0;
};

template <class S_, int dim_>
class ImplicitSphere : public ImplicitSurface<S_, dim_> {
private:
  S_ radius;
  Vector<S_, dim_> center;

public:
  ImplicitSphere(const Vector<S_, dim_> &c, S_ r)
    : radius(r), center(c) {}

  virtual S_ evaluate(const Vector<S_, dim_> &x) const {
    return (x - center).length() - radius;
  }

  virtual Vector<S_, dim_> gradient(const Vector<S_, dim_> &x) const {
    return (x - center) / (x - center).length();
  }
  
  virtual Mat<S_, dim_, dim_> hessian(const Vector<S_, dim_> &x) const {
    S_ dist = (x - center).length();
    S_ dist3 = dist * dist * dist;
    return Mat<S_, dim_, dim_>::identity() / dist;/*(Mat<S_, dim_, dim_>::identity() * (dist * dist) -
	     Mat<S_, dim_, dim_>::outerProduct(x - center, x - center)) / dist3*/;
  }
};

template <class S_, int dim_>
class ImplicitBiharmonicRBF : public ImplicitSurface<S_, dim_> {
private:
  std::vector<Vector<S_, dim_> > centers;
  std::vector<S_> weights;

public:
  ImplicitBiharmonicRBF(const std::vector<Vector<S_, dim_> > &c,
                        const std::vector<S_> &w) 
    : centers(c), weights(w) {}

  virtual S_ evaluate(const Vector<S_, dim_> &x) const {
    S_ f = (S_)0;
    
    for(unsigned int i = 0; i < weights.size(); i++)
      f += (weights[i] * (x - centers[i]).length());
    
    return f;
  }

  virtual Vector<S_, dim_> gradient(const Vector<S_, dim_> &x) const {
    Vector<S_, dim_> g((S_)0);

    for(unsigned int i = 0; i < weights.size(); i++)
      g += ((x - centers[i]).normalize() * weights[i]);
    
    return g;
  }

  virtual Mat<S_, dim_, dim_> hessian(const Vector<S_, dim_> &x) const {
    Mat<S_, dim_, dim_> h((S_)0);
    Mat<S_, dim_, dim_> id = Mat<S_, dim_, dim_>::identity();
    Vector<S_, dim_> ngrad = gradient(x).normalize();

    for(unsigned int i = 0; i < weights.size(); i++)
      h += ((id - Mat<S_, dim_, dim_>::outerProduct(ngrad, ngrad))
            * weights[i] / (x - centers[i]).length());
    
    return h;
  }
};

template <class S_, int dim_>
class BlentSurface : public DynImplicitSurface<S_, dim_> {
private:
  ImplicitSurface<S_, dim_> *first, *second;
  
public:
  BlentSurface(ImplicitSurface<S_, dim_> *f, ImplicitSurface<S_, dim_> *s)
    : first(f), second(s) {}

  virtual S_ evaluate(const Vector<S_, dim_> &x, S_ t) const {
    return first->evaluate(x) * (1 - t) + second->evaluate(x) * t;
  }
  
  virtual Vector<S_, dim_> gradient(const Vector<S_, dim_> &x, S_ t) const {
    return first->gradient(x) * (1 - t) + second->gradient(x) * t;
  }
  
  virtual Mat<S_, dim_, dim_> hessian(const Vector<S_, dim_> &x, S_ t) const {
    return first->hessian(x) * (1 - t) + second->hessian(x) * t;
  }
  
  virtual S_ timeParDer(const Vector<S_, dim_> &x, S_ t) const {
    return second->evaluate(x) - first->evaluate(x);
  }

  virtual Vector<S_, dim_> timeParDerOfGradient(const Vector<S_, dim_> &x,
                                                S_ t) const {
    return second->gradient(x) - first->gradient(x);
  }
};

template <class S_, int dim_>
class DynBlobbySurface : public DynImplicitSurface<S_, dim_> {
public:
  std::vector<S_> radii;
  std::vector<S_> weights;
  std::vector<Vector<S_, dim_> > centers;
  S_ threshold;
  Vector<S_, dim_> velocity;

  DynBlobbySurface(const std::vector<S_> &r, const std::vector<S_> &w,
                   const std::vector<Vector<S_, dim_> > &c, S_ th)
    : radii(r), weights(w), centers(c), threshold(th) {
    assert(radii.size() == weights.size() && radii.size() == centers.size());
    
    velocity = Vector<S_, dim_>(0.0);
    velocity[0] = 0.5;
  }

  virtual S_ evaluate(const Vector<S_, dim_> &x, S_ t) const {
    S_ res = (S_)0;

    for(unsigned int i = 0; i < centers.size(); i++) {
      Vector<S_, dim_> xt = centers[i] + velocity * t;
      res += weights[i] * kernel((x - xt).length() / radii[i]);
    }

    return res - threshold;
  }

  virtual Vector<S_, dim_> gradient(const Vector<S_, dim_> &x, S_ t) const {
    Vector<S_, dim_> grad(0.0);
    
    for(unsigned int i = 0; i < centers.size(); i++) {
      Vector<S_, dim_> xt = centers[i] + velocity * t;
      S_ r = (x - xt).length() / radii[i];
      S_ fact = (1 - r * r);
      
      grad += ((x - xt) * (fact * fact) *
               (weights[i] / (radii[i] * radii[i])));
    }
    
    return grad * (-6.0);
  }
  
  virtual Mat<S_, dim_, dim_> hessian(const Vector<S_, dim_> &x, S_ t) const {
    S_ coeff = 0;

    for(unsigned int i = 0; i < centers.size(); i++) {
      Vector<S_, dim_> xt = centers[i] + velocity * t;
      S_ r = (x - xt).length() / radii[i];
      S_ fact = (1 - r * r);
      
      coeff += ((weights[i] / (radii[i] * radii[i])) * fact * fact);
    }

    return Mat<S_, dim_, dim_>::identity() * (-6.0) * coeff;
  }
  
  virtual S_ timeParDer(const Vector<S_, dim_> &x, S_ t) const {
    S_ dfdt = 0;
    
    for(unsigned int i = 0; i < centers.size(); i++) {
      Vector<S_, dim_> xt = centers[i] + velocity * t;
      S_ r = (x - xt).length() / radii[i];
      S_ fact = (1 - r * r);
      
      dfdt += ((x - xt) * (fact * fact) *
               (weights[i] / (radii[i] * radii[i]))).dot(velocity);
    }
    
    return (6.0) * dfdt;
  }
  
  virtual Vector<S_, dim_> timeParDerOfGradient(const Vector<S_, dim_> &x,
                                                S_ t) const {
    Vector<S_, dim_> dgfdt(0);

    for(unsigned int i = 0; i < centers.size(); i++) {
      Vector<S_, dim_> xt = centers[i] + velocity * t;
      S_ r = (x - xt).length() / radii[i];
      S_ fact = (1 - r * r);

      dgfdt += velocity * ((weights[i] / (radii[i] * radii[i])) * fact * fact);
    }

    return dgfdt * 6.0;
  }
  
private:
  S_ kernel(S_ r) const {
    S_ factor = 1 - r * r;
    return factor * factor * factor;
  }
};

typedef DynImplicitSurface<float, 2> DynImplicitSurface2f;
typedef DynImplicitSurface<float, 3> DynImplicitSurface3f;
typedef DynImplicitSurface<double, 2> DynImplicitSurface2d;
typedef DynImplicitSurface<double, 3> DynImplicitSurface3d;

typedef ImplicitSphere<float, 2> ImplicitSphere2f;
typedef ImplicitSphere<float, 3> ImplicitSphere3f;
typedef ImplicitSphere<double, 2> ImplicitSphere2d;
typedef ImplicitSphere<double, 3> ImplicitSphere3d;

typedef ImplicitBiharmonicRBF<float, 2> ImplicitBiharmonicRBF2f;
typedef ImplicitBiharmonicRBF<float, 3> ImplicitBiharmonicRBF3f;
typedef ImplicitBiharmonicRBF<double, 2> ImplicitBiharmonicRBF2d;
typedef ImplicitBiharmonicRBF<double, 3> ImplicitBiharmonicRBF3d;

typedef BlentSurface<float, 2> BlentSurface2f;
typedef BlentSurface<float, 3> BlentSurface3f;
typedef BlentSurface<double, 2> BlentSurface2d;
typedef BlentSurface<double, 3> BlentSurface3d;

typedef DynBlobbySurface<float, 2> DynBlobbySurface2f;
typedef DynBlobbySurface<float, 3> DynBlobbySurface3f;
typedef DynBlobbySurface<double, 2> DynBlobbySurface2d;
typedef DynBlobbySurface<double, 3> DynBlobbySurface3d;

template<class S_, int dim_>
Vector<S_, dim_> normal_velocity(DynImplicitSurface<S_, dim_> *f,
                                 const Vector<S_, dim_> &x, S_ t) {
  S_ dfdt = f->timeParDer(x, t);
  Vector<S_, dim_> gradf = f->gradient(x, t);
  Vector<S_, dim_> v = gradf * (-dfdt / gradf.dot(gradf));

  return v;
};

template<class S_, int dim_>
Vector<S_, dim_> tangential_velocity(DynImplicitSurface<S_, dim_> *f,
                                     const Vector<S_, dim_> &x, S_ t) {
  Vector<S_, dim_> norm_gradf = f->gradient(x, t).normalize();
  Mat<S_, dim_, dim_> projf = Mat<S_, dim_, dim_>::projection(norm_gradf);
  Mat<S_, dim_, dim_> hess = f->hessian(x, t);

  return -projf * hess.inverse() * f->timeParDerOfGradient(x, t);
};

template<class S_, int dim_>
Vector<S_, dim_> total_velocity(DynImplicitSurface<S_, dim_> *f,
                                const Vector<S_, dim_> &x, S_ t) {
  return normal_velocity(f, x, t) + tangential_velocity(f, x, t);
};

template <class S_, int dim_>
void advect(DynImplicitSurface<S_, dim_> *f, S_ ti, S_ tf,
	    int nsteps, std::vector<Vector<S_, dim_> > &p) {
  S_ t = ti, dt = (tf - ti) / (S_(nsteps));

  while(t < tf) {
    for(unsigned int i = 0; i < p.size(); i++) {
      Vector<S_, dim_> pos = p[i];
      Vector<S_, dim_> vel = total_velocity(f, pos, t);

      p[i] = pos + vel * dt;
    }
    t += dt;
  }
}

template <class S_, int dim_>
ScalarField<S_, dim_> discrete_field(DynImplicitSurface<S_, dim_> *f, S_ t,
				     const GridConf<S_, dim_> &conf) {
  ScalarField<S_, dim_> sfield(conf);

  for(int i = 0; i < sfield.size(); i++) {
    Vector<S_, dim_> x = sfield.coordinateAt(i);
    sfield.at(i) = f->evaluate(x, t);
  }

  return sfield;
}

template <class S_, int dim_>
ScalarField<S_, dim_> discrete_partial_t(DynImplicitSurface<S_, dim_> *f, S_ t,
					 const GridConf<S_, dim_> &conf) {
  ScalarField<S_, dim_> sfield(conf);

  for(int i = 0; i < sfield.size(); i++) {
    Vector<S_, dim_> x = sfield.coordinateAt(i);
    sfield.at(i) = f->timeParDer(x, t);
  }

  return sfield;
}

template <class S_, int dim_>
VectorField<S_, dim_> discrete_grad(DynImplicitSurface<S_, dim_> *f, S_ t,
				    const GridConf<S_, dim_> &conf) {
  VectorField<S_, dim_> vfield(conf);

  for(int i = 0; i < vfield.size(); i++) {
    Vector<S_, dim_> x = vfield.coordinateAt(i);
    vfield.at(i) = f->gradient(x, t);
  }

  return vfield;
}

template <class S_, int dim_>
VectorField<S_, dim_> discrete_partial_t_grad(DynImplicitSurface<S_, dim_> *f, S_ t,
					      const GridConf<S_, dim_> &conf) {
  VectorField<S_, dim_> vfield(conf);

  for(int i = 0; i < vfield.size(); i++) {
    Vector<S_, dim_> x = vfield.coordinateAt(i);
    vfield.at(i) = f->timeParDerOfGradient(x, t);
  }

  return vfield;
}

template <class S_, int dim_>
VectorField<S_, dim_> discrete_partial_grad(DynImplicitSurface<S_, dim_> *f, S_ t,
					    const GridConf<S_, dim_> &conf, int component) {
  VectorField<S_, dim_> vfield(conf);

  for(int i = 0; i < vfield.size(); i++) {
    Vector<S_, dim_> x = vfield.coordinateAt(i);
    vfield.at(i) = f->hessian(x, t).row(component);
  }

  return vfield;
}

template <class S_, int dim_>
VectorField<S_, dim_> discrete_normal_velocity(DynImplicitSurface<S_, dim_> *f, S_ t,
					       const GridConf<S_, dim_> &conf) {
  VectorField<S_, dim_> vfield(conf);

  for(int i = 0; i < vfield.size(); i++) {
    Vector<S_, dim_> x = vfield.coordinateAt(i);
    vfield.at(i) = normal_velocity(f, x, t);
  }

  return vfield;
}

template <class S_, int dim_>
VectorField<S_, dim_> discrete_tangential_velocity(DynImplicitSurface<S_, dim_> *f, S_ t,
						   const GridConf<S_, dim_> &conf) {
  VectorField<S_, dim_> vfield(conf);

  for(int i = 0; i < vfield.size(); i++) {
    Vector<S_, dim_> x = vfield.coordinateAt(i);
    vfield.at(i) = tangential_velocity(f, x, t);
  }

  return vfield;
}

template <class S_, int dim_>
VectorField<S_, dim_> discrete_total_velocity(DynImplicitSurface<S_, dim_> *f, S_ t,
					      const GridConf<S_, dim_> &conf) {
  VectorField<S_, dim_> vfield(conf);

  for(int i = 0; i < vfield.size(); i++) {
    Vector<S_, dim_> x = vfield.coordinateAt(i);
    vfield.at(i) = total_velocity(f, x, t);
  }

  return vfield;
}

#endif
