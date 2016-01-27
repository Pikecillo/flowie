/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __VELOCITY_HPP__
#define __VELOCITY_HPP__

#include <assert.h>

#include "implicit/dynamic_field.hpp"

template <class S_, int dim_>
void surface_normal_velocity(const DynScalarField<S_, dim_> &field,
                             S_ t, S_ dt,
                             VectorField<S_, dim_> &vel) {
  ScalarField<S_, dim_> dfdt(field.getConf());
  VectorField<S_, dim_> gradf(field.getConf());

  dfdt = field.timeDerivative(t, dt);
  gradf = field.grad(t);
  vel = -(gradf * ((dfdt / gradf.dot(gradf))));
}

template <class S_, int dim_>
void surface_tangential_velocity(const DynScalarField<S_, dim_> &field,
                                 S_ t, S_ dt,
                                 VectorField<S_, dim_> &vel) {
  VectorField<S_, dim_> gradf, gradf_dt, n;
  TensorField<S_, dim_> pn, hf;

  // Compute the time derivative of the gradient field
  gradf = field.grad(t);
  gradf_dt = field.timeDerivativeOfGrad(t, dt);
  
  // Compute the normal field
  n = gradf / gradf.magnitude();
  
  // Compute the projection field over the plane perpendicular to normal
  pn = TensorField<S_, dim_>::projection(n);
  
  // Compute the Hessian
  hf = field.hessian(t, dt);

  // Compute tangential velocity
  vel = -(pn * (hf.inverse() * gradf_dt));
}

template <class S_, int dim_>
void surface_total_velocity(const DynScalarField<S_, dim_> &field,
			    S_ t, S_ dt, VectorField<S_, dim_> &vel) {
  VectorField<S_, dim_> norv(field.getConf());
  VectorField<S_, dim_> tanv(field.getConf());

  surface_normal_velocity(field, t, dt, norv);
  surface_tangential_velocity(field, t, dt, tanv);

  vel = norv + tanv;
}

template <class S_, int dim_>
void advect(const VectorField<S_, dim_> &velocity,
            double t, double dt, std::vector<Vector<S_, dim_> > &particles) {
  S_ time = 0;
  Vector<S_, dim_> p, v;

  while(time < t) {
    for(unsigned int i = 0; i < particles.size(); i++) {
      p = particles[i];
      v = velocity.valueAt(p);
      particles[i] = p + v * dt;
    }

    time += dt;
  }
}

template <class S_, int dim_>
void advect(const DynVectorField<S_, dim_> &velocity,
            S_ ti, S_ tf, S_ dt, std::vector<Vector<S_, dim_> > &particles) {
  S_ time = ti;

  while(time <= tf) {
    advect(velocity.at(time, dt), dt, dt, particles);
    time += dt;
  }
}

template <class S_, int dim_>
class DynSurfaceVelocity : public DynVectorField<S_, dim_> {
public:
  DynSurfaceVelocity() {
    dynScalarField = 0;
  }

  DynSurfaceVelocity(const DynSurfaceVelocity<S_, dim_> &other) {
    (*this) = other;
  }

  void operator=(const DynSurfaceVelocity<S_, dim_> &other) {
    this->conf = other.conf();
    dynScalarField = other.dynScalarField;
  }

  void set(const DynScalarField<S_, dim_> *dynScalarF) {
    this->conf = dynScalarF->getConf();
    dynScalarField = dynScalarF;
  }

  virtual VectorField<S_, dim_> at(S_ t, S_ dt) const {
    assert(dynScalarField != 0);

    VectorField<S_, dim_> vel(this->conf);
    surface_total_velocity(*dynScalarField, t, dt, vel);
    return vel;
  }

private:
  GridConf<S_, dim_> conf;
  const DynScalarField<S_, dim_> *dynScalarField;
};

typedef DynSurfaceVelocity<float, 2> DynSurfaceVelocity2f;
typedef DynSurfaceVelocity<float, 3> DynSurfaceVelocity3f;

typedef DynSurfaceVelocity<double, 2> DynSurfaceVelocity2d;
typedef DynSurfaceVelocity<double, 3> DynSurfaceVelocity3d;

#endif
