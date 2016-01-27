/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __FIELD_TYPES__
#define __FIELD_TYPES__

#include <assert.h>

#include <iostream>
#include <vector>

#include "matrix_types.hpp"
#include "vector_types.hpp"

#include "numerical/interpolant.hpp"
#include "utils/reference.hpp"

enum Orientation { XY=0, XZ, YZ };

template <class S_, int dim_>
class GridConf {
private:
  Vector<int, dim_> dim;
  Vector<S_, dim_> lo, hi;

public:
  GridConf() : dim(0), lo((S_)0), hi((S_)1) {}

  GridConf(const Vector<int, dim_> &n,
           const Vector<S_, dim_> &l,
           const Vector<S_, dim_> &h)
    : dim(n), lo(l), hi(h) {
    Vector<int, dim_> two(2); assert(n >= two);
    assert(dim_ >= 2);
    assert(lo <= hi);
  }

  GridConf(const GridConf<S_, dim_> &other) {
    (*this) = other;
  }

  void operator=(const GridConf<S_, dim_> &other) {
    dim = other.dim; lo = other.lo; hi = other.hi;
  }

  bool operator==(const GridConf<S_, dim_> &other) const {
    return (dim == other.dim && lo == other.lo && hi == other.hi);
  }

  Vector<S_, dim_> getLo() const { return lo; }
  Vector<S_, dim_> getHi() const { return hi; }
  Vector<int, dim_> getDim() const { return dim; }
  S_ getLo(int i) const { assert(i >= 0 && i < dim_); return lo[i]; }
  S_ getHi(int i) const { assert(i >= 0 && i < dim_); return hi[i]; }
  int getDim(int i) const { assert(i >= 0 && i < dim_); return dim[i]; }
  int size() const { return dim.volume(); }
  Vector<S_, dim_> sideLength() const { return (hi - lo); };
  Vector<S_, dim_> cellLength() const { return sideLength() / (dim - 1); }
};

template <class S_, class T_, class F_, int dim_>
class BasicField {
protected:
  GridConf<S_, dim_> conf;
  std::vector<T_> cells;

protected:
  void shallowCopy(const F_ &other) {
    conf = other.conf;
    cells = other.cells;
  }

  /*void deepCopy(const F_ &other) {
    conf = other.conf;
    cells.clone(other.cells, conf.size());
    }*/

public:
  Vector<S_, dim_> getLo() const { return conf.getLo(); }
  Vector<S_, dim_> getHi() const { return conf.getHi(); }
  GridConf<S_, dim_> getConf() const { return conf; }
  Vector<int, dim_> getDim() const { return conf.getDim(); }
  int size() const { return conf.size(); }
  Vector<S_, dim_> sideLength() const { return conf.sideLength(); }
  Vector<S_, dim_> cellLength() const { return conf.cellLength(); }

  void clone(const F_ &other) {
    deepCopy(other);
  }

  Vector<int, dim_> stride() const {
    Vector<int, dim_> s;
    s[dim_ - 1] = 1;
    for(int i = dim_ - 2; i >= 0; i--)
      s[i] = conf.getDim(i + 1) * s[i + 1];
    return s;
  }

  Vector<int, dim_> indexAt(int i) const {
    assert(i >= 0 && i < conf.size());
    int count = i;
    Vector<int, dim_> index, s = stride();
    for(int j = 0; j < dim_; j++) {
      index[j] = count / s[j];
      count -= index[j] * s[j];
    }
    return index;
  }

  int indexAt(const Vector<int, dim_> &index) const {
    int i = 0;
    Vector<int, dim_> ixs = index * stride(); 
    for(int j = 0; j < dim_; j++) 
      i += ixs[j];
    return i;
  }

  T_ &at(int i) {
    assert(i >= 0 && i < conf.size());
    return cells[i];
  }

  T_ at(int i) const {
    assert(i >= 0 && i < conf.size());
    return cells[i];
  }

  T_ &at(const Vector<int, dim_> &index) {
    return at(indexAt(index));
  }
  
  T_ at(const Vector<int, dim_> &index) const {
    return at(indexAt(index));
  }

  Vector<S_, dim_> coordinateAt(int i) const {
    return coordinateAt(indexAt(i));
  }

  Vector<S_, dim_> coordinateAt(const Vector<int, dim_> &index) const {
    Vector<int, dim_> zero(0);
    assert(index >= zero && index < getDim());
    return conf.getLo() + index * cellLength();
  }
  
  Vector<S_, dim_> valueAt(const Vector<S_, dim_> &x) const {
    assert(false);
  }

  F_ operator-() const { 
    F_ result(*((F_ *)this));
    for(int i = 0; i < size(); i++)
      result.at(i) = -at(i);
    return result;
  }

  F_ &operator+=(const F_ &other) {
    assert(conf == other.conf);
    for(int i = 0; i < size(); i++)
      at(i) += other.at(i);
    return *((F_ *)this);
  }
  
  F_ &operator-=(const F_ &other) {
    assert(conf == other.conf);
    for(int i = 0; i < size(); i++)
      at(i) -= other.at(i);
    return *((F_ *)this);
  }

  F_ &operator*=(const S_ &scalar) {
    for(int i = 0; i < size(); i++)
      at(i) *= scalar;
    return *((F_ *)this);
  }

  F_ &operator*=(const F_ &other) {
    for(int i = 0; i < size(); i++)
      at(i) *= other.at(i);
    return *((F_ *)this);
  }

  F_ &operator/=(const S_ &scalar) {
    for(int i = 0; i < size(); i++)
      at(i) /= scalar;
    return *((F_ *)this);
  }

  F_ &operator/=(const F_ &other) {
    for(int i = 0; i < size(); i++)
      at(i) /= other.at(i);
    return *((F_ *)this);
  }

  F_ operator+(const F_ &other) const {
    assert(conf == other.conf);
    F_ result(*((F_ *)this));
    result += other;
    return result;
  }

  F_ operator-(const F_ &other) const {
    assert(conf == other.conf);
    F_ result(*((F_ *)this));
    result -= other;
    return result;
  }

  F_ operator*(const S_ &scalar) const {
    F_ result(*((F_ *)this));
    result *= scalar;
    return result;
  }

  F_ operator*(const F_ &other) const {
    assert(conf == other.conf);

    F_ result(*((F_ *)this));
    result *= other;
    return result;
  }

  F_ operator/(const S_ &scalar) const {
    F_ result(*((F_ *)this));
    result /= scalar;
    return result;
  }

  F_ operator/(const F_ &other) const {
    assert(conf == other.conf);
    F_ result(*((F_ *)this));
    result /= other;
    return result;
  }

  F_ linearBlend(const F_ &other, float t) const {
    assert(conf == other.conf);
    F_ blent_field(*((F_ *)this));
    for(int i = 0; i < size(); i++)
      blent_field.at(i) = ((1 - t) * at(i) + t * other.at(i));
    return blent_field;
  }
};

/*template <class S_, class T_, class F_>
T_ BasicField<S_, T_, F_, 2>::valueAt(const Vector<S_, 2> &x) const {
  // Check that value is inside the represented space
  assert(x >= getLo() && p <= getHi());
  
  Vector<S_, 2> length = cellLength();
  Vector<int, 2> index = (x[0] / length.x(), x[1] / length.y());
  Vector<S_, 2> lo, hi;
  
  // Adjust if value is on hi boundary of space
  if(i == this->conf.numCells[0] - 1) i--;
  if(j == this->conf.numCells[1] - 1) j--;
  
  lo = this->coordinateAt(i, j);
  hi = this->coordinateAt(i + 1, j + 1);
  
  Vector<S_, 2> x1(lo.x(), hi.x()), x2(lo.y(), hi.y());
  T_ y_[] = { at(i, j), at(i + 1, j), at(i, j + 1), at(i + 1, j + 1) };
  
  return bilinear_interp(Vector<S_, 2>(x, y), x1, x2, y_);
}

template <class S_, class T_, class F_>
T_ BasicField<S_, T_, F_, 3>::valueAt() const {
  Vector<S_, 3> p(x, y, z);
  
  // Check that value is inside the represented space
  assert(p >= this->conf.gridLo && p <= this->conf.gridHi);
  
  Vector<S_, 3> length = this->sideLength();
  
  int i = (x / length.x()), j = (y / length.y()), k = (z / length.z());
  Vector<S_, 3> lo, hi;
  
  // Adjust if value is on hi boundary of space
  if(i == this->conf.numCells[0] - 1) i--;
  if(j == this->conf.numCells[1] - 1) j--;
  
  lo = coordinateAt(i, j, k);
  hi = coordinateAt(i + 1, j + 1, k + 1);
  
  Vector<S_, 2> x1(lo.x(), hi.x()), x2(lo.y(), hi.y()), x3(lo.z(), hi.z());
  T_ y_[] = { at(i, j, k), at(i + 1, j, k),
	      at(i, j + 1, k), at(i + 1, j + 1, k),
	      at(i, j, k + 1), at(i + 1, j, k + 1),
	      at(i, j + 1, k + 1), at(i + 1, j + 1, k + 1)};
  
  return trilinear_interp(Vector<S_, 3>(x, y, z), x1, x2, x3, y_);
  }*/

template <class S_, int dim_>
class ScalarField : public BasicField<S_, S_, ScalarField<S_, dim_>, dim_> {
public:
  ScalarField() {
    assert(dim_ >= 2);
  }

  ScalarField(const GridConf<S_, dim_> &c) {
    assert(c.size());
    this->conf = c;
    this->cells = std::vector<S_>(c.size());
  }
  
  ScalarField(const ScalarField<S_, dim_> &other) {
    (*this) = other;
  }
  
  void operator=(const ScalarField<S_, dim_> &other) {
    this->shallowCopy(other);
  }
};

template <class S_, int dim_>
class VectorField : public BasicField<S_, Vector<S_, dim_>,
				      VectorField<S_, dim_>, dim_> {
public:
  VectorField() {
    assert(dim_ >= 2);
  }

  VectorField(const GridConf<S_, dim_> &c) {
    this->conf = c;
    this->cells = std::vector<Vector<S_, dim_> >(c.size());
  }
  
  VectorField(const VectorField<S_, dim_> &other) {
    (*this) = other;
  }
  
  void operator=(const VectorField<S_, dim_> &other) {
    this->shallowCopy(other);
  }

  ScalarField<S_, dim_> getComponent(int c) const {
    assert(c >= 0 && c < dim_);
    ScalarField<S_, dim_> component(this->conf);
    for(int i = 0; i < this->size(); i++)
      component.at(i) = this->at(i)[c];
    return component;
  }

  void setComponent(int x, const ScalarField<S_, dim_> &sfield) {
    assert(this->conf == sfield.getConf());
    VectorField<S_, dim_> vfield(this->conf);
    for(int i = 0; i < sfield.size(); i++)
      this->at(i)[x] = sfield.at(i);
  }

  ScalarField<S_, dim_> dot(const VectorField<S_, dim_> &other) const {
    assert(this->conf == other.conf);
    ScalarField<S_, dim_> result(this->conf);
    for(int i = 0; i < this->size(); i++)
      result.at(i) = this->at(i).dot(other.at(i));
    return result;
  }

  ScalarField<S_, dim_> magnitude() const {
    ScalarField<S_, dim_> result(this->conf);
    for(int i = 0; i < this->size(); i++)
      result.at(i) = this->at(i).length();
    return result;
  }
};

template <class S_, int dim_>
class TensorField : public BasicField<S_, Mat<S_, dim_, dim_>,
				      TensorField<S_, dim_>, dim_> {
public:
  TensorField() {
    assert(dim_ >= 2);
  }

  TensorField(const GridConf<S_, dim_> &c) {
    this->conf = c;
    this->cells = std::vector<Mat<S_, dim_, dim_> >(c.size());
  }
  
  TensorField(const TensorField<S_, dim_> &other) {
    (*this) = other;
  }
  
  void operator=(const TensorField<S_, dim_> &other) {
    this->shallowCopy(other);
  }

  VectorField<S_, dim_> getCol(int col) const {
    VectorField<S_, dim_> vfield(this->getConf());
    for(int k = 0; k < this->size(); k++)
      for(int i = 0; i < dim_; i++)
        vfield.at(k)[i] = this->at(k)[i][col];
    return vfield;
  }

  void setCol(int col, const VectorField<S_, dim_> &vfield) {
    assert(this->conf == vfield.getConf());
    for(int k = 0; k < this->size(); k++)
      for(int i = 0; i < dim_; i++)
        this->at(k)[i][col] = vfield.at(k)[i];
  }

  TensorField<S_, dim_> inverse() const {
    TensorField<S_, dim_> result(this->conf);
    for(int i = 0; i < this->size(); i++)
      result.at(i) = this->at(i).inverse();
    return result;
  }

  static TensorField<S_, dim_> identity(const GridConf<S_, dim_> &c) {
    TensorField<S_, dim_> result(c);
    for(int i = 0; i < c.size(); i++)
      result.at(i) = Mat<S_, dim_, dim_>::identity();
    return result;
  }

  static TensorField<S_, dim_> projection(const VectorField<S_, dim_> &n) {
    TensorField<S_, dim_> outer(n.getConf());

    for(int k = 0; k < outer.size(); k++) {
      Mat<S_, dim_, dim_> m;
      Vector<S_, dim_> normal = n.at(k);
      
      for(int i = 0; i < dim_; i++)
        for(int j = 0; j < dim_; j++) {
          m[i][j] = normal[i] * normal[j];
        }
      
      outer.at(k) = m;
    }
    
    return TensorField<S_, dim_>::identity(n.getConf()) - outer;
  }
};

template<class S_, int dim_>
VectorField<S_, dim_> operator*(const VectorField<S_, dim_> &field1,
				const ScalarField<S_, dim_> &field2) {
  assert(field1.getConf() == field2.getConf());
  VectorField<S_, dim_> result(field1.getConf());
  for(int i = 0; i < field1.size(); i++)
    result.at(i) = field1.at(i) * field2.at(i);
  return result;
}

template<class S_, int dim_>
VectorField<S_, dim_> operator*(const TensorField<S_, dim_> &tfield,
				const VectorField<S_, dim_> &vfield) {
  assert(tfield.getConf() == vfield.getConf());
  VectorField<S_, dim_> result(tfield.getConf());
  for(int i = 0; i < tfield.size(); i++)
    result.at(i) = tfield.at(i) * vfield.at(i);
  return result;
}

template<class S_, int dim_>
ScalarField<S_, dim_> operator/(const ScalarField<S_, dim_> &field1,
				const ScalarField<S_, dim_> &field2) {
  assert(field1.getConf() == field2.getConf());
  ScalarField<S_, dim_> result(field1.getConf());
  for(int i = 0; i < field1.size(); i++)
    result.at(i) = field1.at(i) / field2.at(i);
  return result;
}

template<class S_, int dim_>
VectorField<S_, dim_> operator/(const VectorField<S_, dim_> &field1,
				const ScalarField<S_, dim_> &field2) {
  assert(field1.getConf() == field2.getConf());
  VectorField<S_, dim_> result(field1.getConf());
  for(int i = 0; i < field1.size(); i++)
    result.at(i) = field1.at(i) / field2.at(i);
  return result;
}

template <class S_>
void get_slice(const ScalarField<S_, 3> &field3d,
	       ScalarField<S_, 2> &field2d,
	       int orientation, S_ free_coord) {
  GridConf<S_, 3> conf3d = field3d.getConf();
  GridConf<S_, 2> conf2d;

  Vector<int, 3> nc = conf3d.numCells;
  Vector<S_, 3> lo = conf3d.gridLo;
  Vector<S_, 3> hi = conf3d.gridHi;
  if(orientation == XY) {
    conf2d = GridConf<S_, 2>(Vector2i(nc.x(), nc.y()),
			     Vector<S_, 2>(lo.x(), lo.y()),
			     Vector<S_, 2>(hi.x(), hi.y()));
  } else if(orientation == XZ) {
    conf2d = GridConf<S_, 2>(Vector2i(nc.x(), nc.z()),
			     Vector<S_, 2>(lo.x(), lo.z()),
			     Vector<S_, 2>(hi.x(), hi.z()));
  } else {
    conf2d = GridConf<S_, 2>(Vector2i(nc.y(), nc.z()),
			     Vector<S_, 2>(lo.y(), lo.z()),
			     Vector<S_, 2>(hi.y(), hi.z()));
  }

  field2d = ScalarField<S_, 2>(conf2d);

  for(int i = 0; i < conf2d.size(); i++) {
      Vector<S_, 3> p;
      Vector<S_, 2> coord;

      coord = field2d.coordinateAt(i);

      if(orientation == XY) {
	p = Vector<S_, 3>(coord.x(), coord.y(), free_coord);
      } else if(orientation == XZ) {
	p = Vector<S_, 3>(coord.x(), free_coord, coord.y());
      } else {
	p = Vector<S_, 3>(free_coord, coord.x(), coord.y());
      }

      field2d.at(i) = field3d.valueAt(p);
    }
}

typedef GridConf<float, 2> GrifConf2f;
typedef GridConf<float, 3> GrifConf3f;

typedef GridConf<double, 2> GrifConf2d;
typedef GridConf<double, 3> GrifConf3d;

typedef ScalarField<float, 2> ScalarField2f;
typedef ScalarField<float, 3> ScalarField3f;

typedef ScalarField<double, 2> ScalarField2d;
typedef ScalarField<double, 3> ScalarField3d;

typedef VectorField<float, 2> VectorField2f;
typedef VectorField<float, 3> VectorField3f;

typedef VectorField<double, 2> VectorField2d;
typedef VectorField<double, 3> VectorField3d;

typedef TensorField<float, 2> TensorField2f;
typedef TensorField<float, 3> TensorField3f;

typedef TensorField<double, 2> TensorField2d;
typedef TensorField<double, 3> TensorField3d;

#endif
