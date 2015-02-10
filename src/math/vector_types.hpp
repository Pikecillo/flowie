/**
 * Copyright (C) 2013. Mario Rincon-Nigro.
 *
 * This file is a part of Flowie.
 *
 * Flowie is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Flowie is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Flowie.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __VECTOR_TYPES_HPP__
#define __VECTOR_TYPES_HPP__

#include <assert.h>
#include <math.h>

namespace flowie
{

template <class S, template <typename, int> class V, int Size>
class BasicVector {
protected:
    BasicVector() {
	assert(Size >= 1);
    }
    
    explicit BasicVector(S s) {
	assert(Size >= 1);
	for(int i = 0; i < Size; i++)
	    this->elements[i] = s;
    }
    
public:
    int size() const {
	return Size;
    }
    
    S at(int i) const {
	assert(i >= 0 && i < Size);
	return elements[i];
    }
    
    S &at(int i) {
	assert(i >= 0 && i < Size);
	return elements[i];
    }
    
    S operator[](int i) const {
	assert(i >= 0 && i < Size);
	return elements[i];
    }
    
    S &operator[](int i) {
	assert(i >= 0 && i < Size);
	return elements[i];
    }
    
    bool operator==(const V<S, Size> &other) const {
	for(int i = 0; i < Size; i++) {
	    if(elements[i] != other.elements[i])
		return false;
	}
	
	return true;
    }
    
    bool operator==(const S &scalar) const {
	for(int i = 0; i < Size; i++) {
	    if(elements[i] != scalar)
		return false;
	}
	
	return true;
    }
    
    bool operator!=(const V<S, Size> &other) const {
	return !((*this) == other);
    }
    
    bool operator!=(const S &scalar) const {
	return !((*this) == scalar);
    }
    
    V<S, Size> operator-() const {
	V<S, Size> result(*((V<S, Size> *)this));
	
	for(int i = 0; i < Size; i++)
	    result[i] = -elements[i];
	
	return result;
    }
    
    V<S, Size> &operator+=(const S &scalar) {
	for(int i = 0; i < Size; i++)
	    elements[i] += scalar;
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator+=(const V<S, Size> &other) {
	for(int i = 0; i < Size; i++)
	    elements[i] += other.elements[i];
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator-=(const V<S, Size> &other) {
	for(int i = 0; i < Size; i++)
	    elements[i] -= other.elements[i];
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator-=(const S &scalar) {
	for(int i = 0; i < Size; i++)
	    elements[i] -= scalar;
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator*=(const S &scalar) {
	for(int i = 0; i < Size; i++)
	    elements[i] *= scalar;
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator*=(const V<S, Size> &other) {
	for(int i = 0; i < Size; i++)
	    elements[i] *= other.elements[i];
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator/=(const S &scalar) {
	for(int i = 0; i < Size; i++)
	    elements[i] /= scalar;
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> &operator/=(const V<S, Size> &other) {
	for(int i = 0; i < Size; i++)
	    elements[i] /= other.elements[i];
	
	return *((V<S, Size> *)this);
    }
    
    V<S, Size> operator+(const S &scalar) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result += scalar;
	return result;
    }
    
    V<S, Size> operator+(const V<S, Size> &other) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result += other;
	return result;
    }
    
    V<S, Size> operator-(const S &scalar) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result -= scalar;
	return result;
    }
    
    V<S, Size> operator-(const V<S, Size> &other) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result -= other;
	return result;
    }
    
    V<S, Size> operator*(const S &scalar) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result *= scalar;
	return result;
    }
    
    V<S, Size> operator*(const V<S, Size> &other) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result *= other;
	return result;
    }
    
    V<S, Size> operator/(const S &scalar) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result /= scalar;
	return result;
    }
    
    V<S, Size> operator/(const V<S, Size> &other) const {
	V<S, Size> result(*((V<S, Size> *)this));
	result /= other;
	return result;
    }
    
    bool operator<=(const V<S, Size> &other) const {
	for(int i = 0 ; i < Size; i++)
	    if(elements[i] > other.elements[i])
		return false;
	return true;
    }
    
    bool operator>=(const V<S, Size> &other) const {
	for(int i = 0 ; i < Size; i++)
	    if(elements[i] < other.elements[i])
		return false;
	return true;
    }
    
    bool operator<(const V<S, Size> &other) const {
	for(int i = 0 ; i < Size; i++)
	    if(elements[i] >= other.elements[i])
		return false;
	return true;
    }
    
    bool operator>(const V<S, Size> &other) const {
	for(int i = 0 ; i < Size; i++)
	    if(elements[i] <= other.elements[i])
		return false;
	return true;
    }
    
    S dot(const V<S, Size> &other) const {
	S sum = (S)0;
	for(int i = 0; i < Size; i++)
	    sum += elements[i] * other.elements[i];
	return sum;
    }
    
    S length() const {
	V<S, Size> v(*(V<S, Size> *)this);
	return sqrt(v.dot(v));
    }
    
    V<S, Size> normalize() const {
	V<S, Size> unit(*(V<S, Size> *)this);
	unit /= unit.length();
	return unit;
    }
    
    /*
     * Projection of other onto this vector
     */
    V<S, Size> project(const V<S, Size> &other) const {
	V<S, Size> proj(*(V<S, Size> *)this);
	proj = (*this) * (this->dot(other) / this->dot(this));
	return proj;
    }
    
    S volume() const {
	S prod = (S)1;
	for(int i = 0; i < Size; i++)
	    prod *= elements[i];
	return prod;
    }
    
    V<S, Size> max(const V<S, Size> &other) const {
	V<S, Size> result;
	
	for(int i = 0; i < Size; i++) {
	    result[i] = (elements[i] > other.elements[i] ?
			 elements[i] : other.elements[i]);
	}
	
	return result;
    }
    
    V<S, Size> min(const V<S, Size> &other) const {
	V<S, Size> result;
	
	for(int i = 0; i < Size; i++) {
	    result[i] = (elements[i] < other.elements[i] ?
			 elements[i] : other.elements[i]);
	}
	
	return result;
    }
    
    V<S, Size + 1> homogeneous() const {
	V<S, Size + 1> h;
	
	for(int i = 0; i < Size; i++)
	    h[i] = elements[i];
	
	h[Size] = (S)1;
	
	return h;
    }
    
    V<S, Size - 1> cartesian() const {
	V<S, Size - 1> c;
	
	for(int i = 0; i < Size - 1; i++)
	    c[i] = elements[i] / elements[Size - 1];
	
	return c;
    }
    
protected:
    S elements[Size];
};

template <class S, int Size>
class Vector : public BasicVector<S, Vector, Size> {
public:
    Vector() : BasicVector<S, Vector, Size>() {}
    
    explicit Vector(S s)
    : BasicVector<S, Vector, Size>(s) {}
    
    Vector(const Vector<S, Size> &other)
	: BasicVector<S, Vector, Size>() {
	(*this) = other;
    }
    
    void operator=(const Vector<S, Size> &other) {
	for(int i = 0; i < Size; i++)
	    this->elements[i] = other.elements[i];
    }
};

template <class S>
class Vector<S, 2> : public BasicVector<S, Vector, 2> {
public:
    Vector() : BasicVector<S, Vector, 2>() {}
    
    explicit Vector(S s)
    : BasicVector<S, Vector, 2>(s) {}
    
    Vector(const S &e0, const S &e1)
	: BasicVector<S, Vector, 2>() {
	this->elements[0] = e0;
	this->elements[1] = e1;
    }
    
    Vector(const Vector<S, 2> &other)
	: BasicVector<S, Vector, 2>() {
	(*this) = other;
    }
    
    void operator=(const Vector<S, 2> &other) {
	this->elements[0] = other.elements[0];
	this->elements[1] = other.elements[1];
    }
    
    const S &x() const { return this->elements[0]; }
    const S &y() const { return this->elements[1]; }
    S &x() { return this->elements[0]; }
    S &y() { return this->elements[1]; }
};

template <class S>
class Vector<S, 3> : public BasicVector<S, Vector, 3> {
public:
    Vector() : BasicVector<S, Vector, 3>() {}
    
    explicit Vector(S s)
    : BasicVector<S, Vector, 3>(s) {}
    
  Vector(const S &e0, const S &e1, const S &e2)
      : BasicVector<S, Vector, 3>() {
      this->elements[0] = e0;
      this->elements[1] = e1;
      this->elements[2] = e2;
  }
    
    Vector(const Vector<S, 3> &other)
	: BasicVector<S, Vector, 3>() {
	(*this) = other;
    }
    
    void operator=(const Vector<S, 3> &other) {
	this->elements[0] = other.elements[0];
	this->elements[1] = other.elements[1];
	this->elements[2] = other.elements[2];
    }
    
    const S &x() const { return this->elements[0]; }
    const S &y() const { return this->elements[1]; }
    const S &z() const { return this->elements[2]; }
    S &x() { return this->elements[0]; }
    S &y() { return this->elements[1]; }
    S &z() { return this->elements[2]; }
    
    Vector<S, 3> cross(const Vector<S, 3> &other) const {
	Vector<S, 3> res;
	
	res.x() = y() * other.z() - z() * other.y();
	res.y() = z() * other.x() - x() * other.z();
	res.z() = x() * other.y() - y() * other.x();
	
	return res;
    }
};

template <class S>
class Vector<S, 4> : public BasicVector<S, Vector, 4> {
public:
    Vector() : BasicVector<S, Vector, 4>() {}
    
    explicit Vector(S s)
    : BasicVector<S, Vector, 4>(s) {}
    
    Vector(const S &e0, const S &e1, const S &e2, const S &e3)
	: BasicVector<S, Vector, 4>() {
	this->elements[0] = e0;
	this->elements[1] = e1;
	this->elements[2] = e2;
	this->elements[3] = e3;
    }
    
    Vector(const Vector<S, 4> &other)
	: BasicVector<S, Vector, 4>() {
	(*this) = other;
    }
    
    void operator=(const Vector<S, 4> &other) {
	this->elements[0] = other.elements[0];
	this->elements[1] = other.elements[1];
	this->elements[2] = other.elements[2];
	this->elements[3] = other.elements[3];
    }
    
    const S &x() const { return this->elements[0]; }
    const S &y() const { return this->elements[1]; }
    const S &z() const { return this->elements[2]; }
    const S &w() const { return this->elements[3]; }
    S &x() { return this->elements[0]; }
    S &y() { return this->elements[1]; }
    S &z() { return this->elements[2]; }
    S &w() { return this->elements[3]; }
};

#define SCALAR_OP_VECTOR_TEMPLATE(OP, DIM)				\
    template <class S>							\
    inline Vector<S, DIM> operator OP (S s, const Vector<S, DIM> &v) {	\
	Vector<S, DIM> res;						\
	for(int i = 0; i < DIM; i++)					\
	    res[i] = s OP v[i];						\
	return res;							\
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
    template <class A, class B, class C, int Size>	       \
    inline C operator OP (const A &a, const B &b) {	       \
	C res;						       \
	for(int i = 0; i < Size; i++)			       \
	    res[i] = a[i] OP b[i];			       \
	return res;					       \
    }							       \
    
#define OPERATOR_DEFINITION(RET_TYPE, TYPE_A, TYPE_B, OP, DIM)		\
    inline Vector<RET_TYPE, DIM> operator OP (const Vector<TYPE_A, DIM> &a, \
					      const Vector<TYPE_B, DIM> &b) { \
	return operator OP <Vector<TYPE_A, DIM>, Vector<TYPE_B, DIM>,	\
			    Vector<RET_TYPE, DIM>, DIM>(a, b);		\
    }									\
    
#define ALL_PAIRS_MIXED_TYPE_OPERATOR_DEFINITION(OP)			\
    OPERATOR_DEFINITION(float, float, int, OP, 2)			\
    OPERATOR_DEFINITION(float, int, float, OP, 2)			\
    OPERATOR_DEFINITION(double, double, int, OP, 2)			\
    OPERATOR_DEFINITION(double, int, double, OP, 2)			\
    OPERATOR_DEFINITION(double, double, float, OP, 2)			\
    OPERATOR_DEFINITION(double, float, double, OP, 2)			\
    OPERATOR_DEFINITION(float, float, int, OP, 3)			\
    OPERATOR_DEFINITION(float, int, float, OP, 3)			\
    OPERATOR_DEFINITION(double, double, int, OP, 3)			\
    OPERATOR_DEFINITION(double, int, double, OP, 3)			\
    OPERATOR_DEFINITION(double, double, float, OP, 3)			\
    OPERATOR_DEFINITION(double, float, double, OP, 3)			\
    OPERATOR_DEFINITION(float, float, int, OP, 4)			\
    OPERATOR_DEFINITION(float, int, float, OP, 4)			\
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

}

#endif
