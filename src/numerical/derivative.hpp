#ifndef __DERIVATIVE_HPP__
#define __DERIVATIVE_HPP__

#include <assert.h>

#include "field_types.hpp"

template <class T_, int dim_>
void partial_derivative(int x, const ScalarField<T_, dim_> &field,
                        ScalarField<T_, dim_> &derivative) {
  assert(x >= 0 && x < dim_);
  assert(field.getConf() == derivative.getConf());

  Vector<T_, dim_> h = field.cellLength();
  Vector<int, dim_> s = field.stride(), size = field.getDim();

  for(int i = 0; i < derivative.size(); i++) {
    // The derivatives are the average of the forward
    // and backward differences
    
    Vector<int, dim_> index = field.indexAt(i);
    T_ fbck, ffwd;

    // Forward/backward values needed at boundaries are extrapolated
    if(index[x] == 0) {
      fbck = 2 * field.at(i) - field.at(i + s[x]);
    } else {
      fbck = field.at(i - s[x]);
    }
    if(index[x] == size[x] - 1) {
      ffwd = 2 * field.at(i) - field.at(i - s[x]);
    } else {
      ffwd = field.at(i + s[x]);
    }

    derivative.at(i) = (ffwd - fbck) / (2 * h[x]);
  }
}

template <class T_, int dim_>
void gradient(const ScalarField<T_, dim_> &sfield,
              VectorField<T_, dim_> &vfield) {
  assert(sfield.getConf() == vfield.getConf());

  ScalarField<T_, dim_> derivative(sfield.getConf());
  
  // For each dimension
  for(int x = 0; x < dim_; x++) {
    // Compute the partial derivative
    partial_derivative(x, sfield, derivative);
    // And set the corresponding component in the gradient field
    vfield.setComponent(x, derivative);
  }
}

template<class S_, int dim_>
void hessian(const ScalarField<S_, dim_> &f,
	     TensorField<S_, dim_> &h) {
  assert(f.getConf() == h.getConf());

  VectorField<S_, dim_> grad_f(f.getConf());
  VectorField<S_, dim_> grad_dfdx(f.getConf());

  gradient(f, grad_f);

  for(int i = 0; i < dim_; i++) {
    VectorField<S_, dim_> grad_dfdx(f.getConf());
    gradient(grad_f.getComponent(i), grad_dfdx);
    h.setCol(i, grad_dfdx);
  }
}

#endif
