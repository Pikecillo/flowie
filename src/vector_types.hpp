/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __VECTOR_TYPES_HPP__
#define __VECTOR_TYPES_HPP__

#include <assert.h>
#include <math.h>

template <class S_, template <typename, int> class V_, int size_>
class BasicVector {
protected:
  BasicVector() {
    assert(size_ >= 1);
  }

  explicit BasicVector(S_ s) {
    assert(size_ >= 1);
    for(int i = 0; i < size_; i++)
      this->elements[i] = s;
  }

public:
  int size() const {
    return size_;
  }

  S_ at(int i) const {
    assert(i >= 0 && i < size_);
    return elements[i];
  }
  
  S_ &at(int i) {
    assert(i >= 0 && i < size_);
    return elements[i];
  }

  S_ operator[](int i) const {
    assert(i >= 0 && i < size_);
    return elements[i];
  }

  S_ &operator[](int i) {
    assert(i >= 0 && i < size_);
    return elements[i];
  }

  bool operator==(const V_<S_, size_> &other) const {
    for(int i = 0; i < size_; i++) {
      if(elements[i] != other.elements[i])
        return false;
    }

    return true;
  }

  bool operator==(const S_ &scalar) const {
    for(int i = 0; i < size_; i++) {
      if(elements[i] != scalar)
        return false;
    }

    return true;
  }

  bool operator!=(const V_<S_, size_> &other) const {
    return !((*this) == other);
  }

  bool operator!=(const S_ &scalar) const {
    return !((*this) == scalar);
  }

  V_<S_, size_> operator-() const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    
    for(int i = 0; i < size_; i++)
      result[i] = -elements[i];

    return result;
  }

  V_<S_, size_> &operator+=(const S_ &scalar) {
    for(int i = 0; i < size_; i++)
      elements[i] += scalar;

    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> &operator+=(const V_<S_, size_> &other) {
    for(int i = 0; i < size_; i++)
      elements[i] += other.elements[i];

    return *((V_<S_, size_> *)this);
  }
  
  V_<S_, size_> &operator-=(const V_<S_, size_> &other) {
    for(int i = 0; i < size_; i++)
      elements[i] -= other.elements[i];

    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> &operator-=(const S_ &scalar) {
    for(int i = 0; i < size_; i++)
      elements[i] -= scalar;

    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> &operator*=(const S_ &scalar) {
    for(int i = 0; i < size_; i++)
      elements[i] *= scalar;
    
    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> &operator*=(const V_<S_, size_> &other) {
    for(int i = 0; i < size_; i++)
      elements[i] *= other.elements[i];

    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> &operator/=(const S_ &scalar) {
    for(int i = 0; i < size_; i++)
      elements[i] /= scalar;
    
    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> &operator/=(const V_<S_, size_> &other) {
    for(int i = 0; i < size_; i++)
      elements[i] /= other.elements[i];

    return *((V_<S_, size_> *)this);
  }

  V_<S_, size_> operator+(const S_ &scalar) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result += scalar;
    return result;
  }

  V_<S_, size_> operator+(const V_<S_, size_> &other) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result += other;
    return result;
  }

  V_<S_, size_> operator-(const S_ &scalar) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result -= scalar;
    return result;
  }

  V_<S_, size_> operator-(const V_<S_, size_> &other) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result -= other;
    return result;
  }

  V_<S_, size_> operator*(const S_ &scalar) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result *= scalar;
    return result;
  }

  V_<S_, size_> operator*(const V_<S_, size_> &other) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result *= other;
    return result;
  }

  V_<S_, size_> operator/(const S_ &scalar) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result /= scalar;
    return result;
  }

  V_<S_, size_> operator/(const V_<S_, size_> &other) const {
    V_<S_, size_> result(*((V_<S_, size_> *)this));
    result /= other;
    return result;
  }

  bool operator<=(const V_<S_, size_> &other) const {
    for(int i = 0 ; i < size_; i++)
      if(elements[i] > other.elements[i])
	return false;
    return true;
  }

  bool operator>=(const V_<S_, size_> &other) const {
    for(int i = 0 ; i < size_; i++)
      if(elements[i] < other.elements[i])
	return false;
    return true;
  }

   bool operator<(const V_<S_, size_> &other) const {
    for(int i = 0 ; i < size_; i++)
      if(elements[i] >= other.elements[i])
	return false;
    return true;
  }

  bool operator>(const V_<S_, size_> &other) const {
    for(int i = 0 ; i < size_; i++)
      if(elements[i] <= other.elements[i])
	return false;
    return true;
  }

  S_ dot(const V_<S_, size_> &other) const {
    S_ sum = (S_)0;
    for(int i = 0; i < size_; i++)
      sum += elements[i] * other.elements[i];
    return sum;
  }

  S_ length() const {
    V_<S_, size_> v(*(V_<S_, size_> *)this);
    return sqrt(v.dot(v));
  }

  V_<S_, size_> normalize() const {
    V_<S_, size_> unit(*(V_<S_, size_> *)this);
    unit /= unit.length();
    return unit;
  }

  /*
   * Projection of other onto this vector
   */
  V_<S_, size_> project(const V_<S_, size_> &other) const {
    V_<S_, size_> proj(*(V_<S_, size_> *)this);
    proj = (*this) * (this->dot(other) / this->dot(this));
    return proj;
  }

  S_ volume() const {
    S_ prod = (S_)1;
    for(int i = 0; i < size_; i++)
      prod *= elements[i];
    return prod;
  }

  V_<S_, size_> max(const V_<S_, size_> &other) const {
    V_<S_, size_> result;
    
    for(int i = 0; i < size_; i++) {
      result[i] = (elements[i] > other.elements[i] ?
		   elements[i] : other.elements[i]);
    }

    return result;
  }

  V_<S_, size_> min(const V_<S_, size_> &other) const {
    V_<S_, size_> result;
    
    for(int i = 0; i < size_; i++) {
      result[i] = (elements[i] < other.elements[i] ?
		   elements[i] : other.elements[i]);
    }

    return result;
  }

  V_<S_, size_ + 1> homogeneous() const {
    V_<S_, size_ + 1> h;

    for(int i = 0; i < size_; i++)
      h[i] = elements[i];

    h[size_] = (S_)1;

    return h;
  }

  V_<S_, size_ - 1> cartesian() const {
    V_<S_, size_ - 1> c;

    for(int i = 0; i < size_ - 1; i++)
      c[i] = elements[i] / elements[size_ - 1];

    return c;
  }

protected:
  S_ elements[size_];
};

template <class S_, int size_>
class Vector : public BasicVector<S_, Vector, size_> {
public:
  Vector() : BasicVector<S_, Vector, size_>() {}

  explicit Vector(S_ s)
    : BasicVector<S_, Vector, size_>(s) {}

  Vector(const Vector<S_, size_> &other)
    : BasicVector<S_, Vector, size_>() {
    (*this) = other;
  }

  void operator=(const Vector<S_, size_> &other) {
    for(int i = 0; i < size_; i++)
      this->elements[i] = other.elements[i];
  }
};

template <class S_>
class Vector<S_, 2> : public BasicVector<S_, Vector, 2> {
public:
  Vector() : BasicVector<S_, Vector, 2>() {}

  explicit Vector(S_ s)
    : BasicVector<S_, Vector, 2>(s) {}

  Vector(const S_ &e0, const S_ &e1)
    : BasicVector<S_, Vector, 2>() {
    this->elements[0] = e0;
    this->elements[1] = e1;
  }

  Vector(const Vector<S_, 2> &other)
    : BasicVector<S_, Vector, 2>() {
    (*this) = other;
  }
  
  void operator=(const Vector<S_, 2> &other) {
    this->elements[0] = other.elements[0];
    this->elements[1] = other.elements[1];
  }

  const S_ &x() const { return this->elements[0]; }
  const S_ &y() const { return this->elements[1]; }
  S_ &x() { return this->elements[0]; }
  S_ &y() { return this->elements[1]; }
};

template <class S_>
class Vector<S_, 3> : public BasicVector<S_, Vector, 3> {
public:
  Vector() : BasicVector<S_, Vector, 3>() {}

  explicit Vector(S_ s)
    : BasicVector<S_, Vector, 3>(s) {}

  Vector(const S_ &e0, const S_ &e1, const S_ &e2)
    : BasicVector<S_, Vector, 3>() {
    this->elements[0] = e0;
    this->elements[1] = e1;
    this->elements[2] = e2;
  }

  Vector(const Vector<S_, 3> &other)
    : BasicVector<S_, Vector, 3>() {
    (*this) = other;
  }
  
  void operator=(const Vector<S_, 3> &other) {
    this->elements[0] = other.elements[0];
    this->elements[1] = other.elements[1];
    this->elements[2] = other.elements[2];
  }

  const S_ &x() const { return this->elements[0]; }
  const S_ &y() const { return this->elements[1]; }
  const S_ &z() const { return this->elements[2]; }
  S_ &x() { return this->elements[0]; }
  S_ &y() { return this->elements[1]; }
  S_ &z() { return this->elements[2]; }

  Vector<S_, 3> cross(const Vector<S_, 3> &other) const {
    Vector<S_, 3> res;
    
    res.x() = y() * other.z() - z() * other.y();
    res.y() = z() * other.x() - x() * other.z();
    res.z() = x() * other.y() - y() * other.x();
    
    return res;
  }
};

template <class S_>
class Vector<S_, 4> : public BasicVector<S_, Vector, 4> {
public:
  Vector() : BasicVector<S_, Vector, 4>() {}

  explicit Vector(S_ s)
    : BasicVector<S_, Vector, 4>(s) {}

  Vector(const S_ &e0, const S_ &e1, const S_ &e2, const S_ &e3)
    : BasicVector<S_, Vector, 4>() {
    this->elements[0] = e0;
    this->elements[1] = e1;
    this->elements[2] = e2;
    this->elements[3] = e3;
  }

  Vector(const Vector<S_, 4> &other)
    : BasicVector<S_, Vector, 4>() {
    (*this) = other;
  }
  
  void operator=(const Vector<S_, 4> &other) {
    this->elements[0] = other.elements[0];
    this->elements[1] = other.elements[1];
    this->elements[2] = other.elements[2];
    this->elements[3] = other.elements[3];
  }

  const S_ &x() const { return this->elements[0]; }
  const S_ &y() const { return this->elements[1]; }
  const S_ &z() const { return this->elements[2]; }
  const S_ &w() const { return this->elements[3]; }
  S_ &x() { return this->elements[0]; }
  S_ &y() { return this->elements[1]; }
  S_ &z() { return this->elements[2]; }
  S_ &w() { return this->elements[3]; }
};

#define SCALAR_OP_VECTOR_TEMPLATE(OP, DIM)				\
  template <class S_>							\
  inline Vector<S_, DIM> operator OP (S_ s, const Vector<S_, DIM> &v) {	\
    Vector<S_, DIM> res;						\
    for(int i = 0; i < DIM; i++)					\
      res[i] = s OP v[i];						\
    return res;								\
  }									\

SCALAR_OP_VECTOR_TEMPLATE(+, 2)
SCALAR_OP_VECTOR_TEMPLATE(+, 3)
SCALAR_OP_VECTOR_TEMPLATE(+, 4)
SCALAR_OP_VECTOR_TEMPLATE(-, 2)
SCALAR_OP_VECTOR_TEMPLATE(-, 3)
SCALAR_OP_VECTOR_TEMPLATE(-, 4)
SCALAR_OP_VECTOR_TEMPLATE(*, 2)
SCALAR_OP_VECTOR_TEMPLATE(*, 3)
SCALAR_OP_VECTOR_TEMPLATE(*, 4)
SCALAR_OP_VECTOR_TEMPLATE(/, 2)
SCALAR_OP_VECTOR_TEMPLATE(/, 3)
SCALAR_OP_VECTOR_TEMPLATE(/, 4)

#define MIXED_TYPE_OPERATOR_TEMPLATE(OP)                       \
  template <class A_, class B_, class C_, int size_>	       \
  inline C_ operator OP (const A_ &a, const B_ &b) {	       \
    C_ res;						       \
    for(int i = 0; i < size_; i++)			       \
      res[i] = a[i] OP b[i];				       \
    return res;						       \
  }							       \
  
#define OPERATOR_DEFINITION(RET_TYPE, TYPE_A, TYPE_B, OP, DIM)		\
  inline Vector<RET_TYPE, DIM> operator OP (const Vector<TYPE_A, DIM> &a, \
					    const Vector<TYPE_B, DIM> &b) { \
    return operator OP <Vector<TYPE_A, DIM>, Vector<TYPE_B, DIM>,	\
			Vector<RET_TYPE, DIM>, DIM>(a, b);		\
  }									\
  
#define ALL_PAIRS_MIXED_TYPE_OPERATOR_DEFINITION(OP)			\
  OPERATOR_DEFINITION(float, float, int, OP, 2)				\
  OPERATOR_DEFINITION(float, int, float, OP, 2)				\
  OPERATOR_DEFINITION(double, double, int, OP, 2)			\
  OPERATOR_DEFINITION(double, int, double, OP, 2)			\
  OPERATOR_DEFINITION(double, double, float, OP, 2)			\
  OPERATOR_DEFINITION(double, float, double, OP, 2)			\
  OPERATOR_DEFINITION(float, float, int, OP, 3)				\
  OPERATOR_DEFINITION(float, int, float, OP, 3)				\
  OPERATOR_DEFINITION(double, double, int, OP, 3)			\
  OPERATOR_DEFINITION(double, int, double, OP, 3)			\
  OPERATOR_DEFINITION(double, double, float, OP, 3)			\
  OPERATOR_DEFINITION(double, float, double, OP, 3)			\
  OPERATOR_DEFINITION(float, float, int, OP, 4)				\
  OPERATOR_DEFINITION(float, int, float, OP, 4)				\
  OPERATOR_DEFINITION(double, double, int, OP, 4)			\
  OPERATOR_DEFINITION(double, int, double, OP, 4)			\
  OPERATOR_DEFINITION(double, double, float, OP, 4)			\
  OPERATOR_DEFINITION(double, float, double, OP, 4)			\
  
MIXED_TYPE_OPERATOR_TEMPLATE(+)
MIXED_TYPE_OPERATOR_TEMPLATE(-)
MIXED_TYPE_OPERATOR_TEMPLATE(/)
MIXED_TYPE_OPERATOR_TEMPLATE(*)

ALL_PAIRS_MIXED_TYPE_OPERATOR_DEFINITION(+)
ALL_PAIRS_MIXED_TYPE_OPERATOR_DEFINITION(-)
ALL_PAIRS_MIXED_TYPE_OPERATOR_DEFINITION(*)
ALL_PAIRS_MIXED_TYPE_OPERATOR_DEFINITION(/)

typedef Vector<int, 2> Vector2i;
typedef Vector<int, 3> Vector3i;
typedef Vector<int, 4> Vector4i;

typedef Vector<float, 2> Vector2f;
typedef Vector<float, 3> Vector3f;
typedef Vector<float, 4> Vector4f;

typedef Vector<double, 2> Vector2d;
typedef Vector<double, 3> Vector3d;
typedef Vector<double, 4> Vector4d;

typedef Vector2i Size2i;
typedef Vector3i Size3i;
typedef Vector4i Size4i;

typedef Vector2f Size2f;
typedef Vector3f Size3f;
typedef Vector4f Size4f;

typedef Vector2d Size2d;
typedef Vector3d Size3d;
typedef Vector4d Size4d;

#endif
