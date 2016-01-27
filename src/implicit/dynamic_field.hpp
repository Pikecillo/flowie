/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __DYNAMIC_FIELD_HPP__
#define __DYNAMIC_FIELD_HPP__

#include <vector>

#include "field_types.hpp"
#include "numerical/derivative.hpp"

template <class S_, int dim_>
class DynScalarField {
public:
  virtual ScalarField<S_, dim_> at(S_ t) const = 0;
  virtual S_ at(const Vector<int, dim_> &index, S_ t) const = 0;
  virtual ScalarField<S_, dim_> timeDerivative(S_ t, S_ dt) const = 0;
  virtual VectorField<S_, dim_> grad(S_ t) const = 0;
  virtual VectorField<S_, dim_> timeDerivativeOfGrad(S_ t, S_ dt) const = 0;
  virtual TensorField<S_, dim_> hessian(S_ t, S_ dt) const = 0;

  Vector<int, dim_> size() const { return conf.numCells; }
  Vector<S_, dim_> getGridLo() const { return conf.gridLo; }
  Vector<S_, dim_> getGridHi() const { return conf.gridHi; }
  GridConf<S_, dim_> getConf() const { return conf; }

protected:
  GridConf<S_, dim_> conf;
};

template <class S_, int dim_>
class DynVectorField {
public:
  virtual VectorField<S_, dim_> at(S_ t, S_ dt) const = 0;

  Vector<int, dim_> size() const { return conf.numCells; }
  Vector<S_, dim_> getGridLo() const { return conf.gridLo; }
  Vector<S_, dim_> getGridHi() const { return conf.gridHi; }
  GridConf<S_, dim_> getConf() const { return conf; }

protected:
  GridConf<S_, dim_> conf;
};

template <class S_, int dim_>
class BlendScalarField : public DynScalarField<S_, dim_> {
public:
  BlendScalarField() {}

  void set(const ScalarField<S_, dim_> *f0,
           const ScalarField<S_, dim_> *f1) {
    assert(f0->getConf() == f1->getConf());
    
    field0 = f0; field1 = f1;
    this->conf = f0->getConf();
    
    gfield0 = VectorField<S_, dim_>(f0->getConf());
    gfield1 = VectorField<S_, dim_>(f1->getConf());
    hfield0 = TensorField<S_, dim_>(f0->getConf());
    hfield1 = TensorField<S_, dim_>(f1->getConf());

    ::gradient(*field0, gfield0);
    ::gradient(*field1, gfield1);
    ::hessian(*field0, hfield0);
    ::hessian(*field1, hfield1);
  }

  virtual ScalarField<S_, dim_> at(S_ t) const {
    return field0->linearBlend(*field1, t);
  }

  virtual S_ at(const Vector<int, dim_> &index, S_ t) const {
    return field0->at(index) * (1 - t) + field1->at(index) * t;
  }

  virtual ScalarField<S_, dim_> timeDerivative(S_ t, S_ dt) const {
    return (*field1) - (*field0);
  }

  virtual VectorField<S_, dim_> grad(S_ t) const {
    return gfield0 * (1 - t) + gfield1 * t;
  }

  virtual VectorField<S_, dim_> timeDerivativeOfGrad(S_ t, S_ dt) const {
    return gfield1 - gfield0;
  }

  virtual TensorField<S_, dim_> hessian(S_ t, S_ dt) const {
    return hfield0 * (1 - t) + hfield1 * t;
  }

private:
  const ScalarField<S_, dim_> *field0, *field1;
  VectorField<S_, dim_> gfield0, gfield1;
  TensorField<S_, dim_> hfield0, hfield1;
};

typedef BlendScalarField<float, 2> BlendScalarField2f;
typedef BlendScalarField<float, 3> BlendScalarField3f;

typedef BlendScalarField<double, 2> BlendScalarField2d;
typedef BlendScalarField<double, 2> BlendScalarField3d;

#endif
