/*
 * University of Houston
 * Mario Rincon-Nigro. April 2013.
 */

#ifndef __MATRIX_TYPES_HPP__
#define __MATRIX_TYPES_HPP__

#include <math.h>

#include <iostream>

#include "vector_types.hpp"

template <template <typename, int, int> class M_, typename S_,
	  int rows_, int cols_>
class BasicMat {
protected:
  void set(S_ s) {
    for(int i = 0; i < rows_; i++)
      elements[i] = Vector<S_, cols_>(s);
  }
  
public:
  const Vector<S_, cols_> operator[](int i) const {
    assert(i >= 0 && i <= rows_);
    return elements[i];
  }
  
  Vector<S_, cols_> &operator[](int i) {
    assert(i >= 0 && i <= rows_);
    return elements[i];
  }
  
  const Vector<S_, cols_> row(int i) const {
    assert(i >= 0 && i <= rows_);
    return elements[i];
  }
  
  Vector<S_, cols_> &row(int i) {
    assert(i >= 0 && i <= rows_);
    return elements[i];
  }
 
  Vector<S_, rows_> col(int j) const {
    assert(j >= 0 && j <= cols_);
    Vector<S_, rows_> c;
    for(int i = 0; i < rows_; i++)
      c[i] = elements[i][j];
    return c;
  }

  M_<S_, rows_, cols_> &operator+=(const M_<S_, rows_, cols_> &other) {
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        elements[i][j] += other.elements[i][j];
    
    return *((M_<S_, rows_, cols_> *)this);
  }
  
  M_<S_, rows_, cols_> &operator-=(const M_<S_, rows_, cols_> &other) {
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        elements[i][j] -= other.elements[i][j];
    
    return *((M_<S_, rows_, cols_> *)this);
  }
  
  M_<S_, rows_, cols_> &operator*=(const S_ &scalar) {
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        elements[i][j] *= scalar;
    
    return *((M_<S_, rows_, cols_> *)this);
  }
  
  template <int rows_o, int cols_o>
  M_<S_, rows_, cols_o> &operator*=(const M_<S_, rows_o, cols_o> &other) {
    M_<S_, rows_, cols_o> result(*((M_<S_, rows_, cols_o> *)this));
    
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_o; j++) {
        S_ dot = 0.0;
	
        for(int k = 0; k < cols_; k++)
          dot += (result.elements[i][k] * other.elements[k][j]);
	
        elements[i][j] = dot;
      }
    
    return *((M_<S_, rows_, cols_o> *)this); 
  }
  
  M_<S_, rows_, cols_> &operator/=(const S_ &scalar) {
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        elements[i][j] /= scalar;
    
    return *((M_<S_, rows_, cols_> *)this);
  }

  M_<S_, rows_, cols_> operator+(const M_<S_, rows_, cols_> &other) const {
    M_<S_, rows_, cols_> result(*((M_<S_, rows_, cols_> *)this));
    result += other;
    return result;
  }

  M_<S_, rows_, cols_> operator-() const {
    M_<S_, rows_, cols_> result(*((M_<S_, rows_, cols_> *)this));
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        result[i][j] = -elements[i][j];
    
    return result;
  }
  
  M_<S_, rows_, cols_> operator-(const M_<S_, rows_, cols_> &other) const {
    M_<S_, rows_, cols_> result(*((M_<S_, rows_, cols_> *)this));
    result -= other;
    return result;
  }

  M_<S_, rows_, cols_> operator*(const S_ &scalar) const {
    M_<S_, rows_, cols_> result(*((M_<S_, rows_, cols_> *)this));
    result *= scalar;
    return result;
  }

  Vector<S_, rows_> operator*(const Vector<S_, cols_> &vector) const {
    Vector<S_, rows_> result;
    for(int i = 0; i < rows_; i++)
      result[i] = elements[i].dot(vector);
    
    return result;
  }

  template <int rows_o, int cols_o>
  M_<S_, rows_, cols_o> operator*
  (const M_<S_, rows_o, cols_o> &other) const {
    M_<S_, rows_, cols_o> result(*((M_<S_, rows_, cols_o> *)this));
    result *= other;
    return result;
  }
  
  M_<S_, rows_, cols_> operator/(const S_ &scalar) const {
    M_<S_, rows_, cols_> result(*((M_<S_, rows_, cols_> *)this));
    result /= scalar;
    return result;
  }

  M_<S_, cols_, rows_> transpose() const {
    M_<S_, cols_, rows_> trans(*((M_<S_, cols_, rows_> *)this));
    
    for(int i = 0; i < cols_; i++)
      for(int j = 0; j < rows_; j++)
        trans[i][j] = this->elements[j][i];
    
    return trans;
  }
  
  M_<S_, rows_ - 1, cols_ - 1> subMatrix(int i, int j) const {
    assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
    
    M_<S_, rows_ - 1, cols_ - 1> m;
    for(int x = 0; x < rows_ - 1; x++)
      for(int y = 0; y < cols_ - 1; y++)
	m[x][y] = this->elements[(x < i ? x : x + 1)][(y < j ? y : y + 1)];
    return m;
  }
  
  static M_<S_, rows_, cols_> identity() {
    M_<S_, rows_, cols_> id;
    
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        id[i][j] = (i == j ? (S_)1 : (S_)0);
    
    return id;
  }
  
  template <int dim_>
  static M_<S_, dim_, dim_> outerProduct(const Vector<S_, dim_> &a,
					 const Vector<S_, dim_> &b) {
    M_<S_, dim_, dim_> outer;
    
    for(int i = 0; i < rows_; i++)
      for(int j = 0; j < cols_; j++)
        outer[i][j] = a[i] * b[j];
    
    return outer;
  }
  
  template <int dim_>
  static M_<S_, dim_, dim_> projection(const Vector<S_, dim_> &v) {
    return (M_<S_, dim_, dim_>::identity() -
	    M_<S_, dim_, dim_>::outerProduct(v, v));
  }
  
protected:
  Vector<S_, cols_> elements[rows_]; 
};

template <typename S_, int rows_, int cols_>
class Mat : public BasicMat<Mat, S_, rows_, cols_> {
public:
  Mat() { assert(rows_ > 0 && cols_ > 0); }
 
  Mat(S_ s) { this->set(s); }

  Mat(const Mat<S_, rows_, cols_> &other) {
    (*this) = other;
  }

  Mat<S_, rows_, cols_> &operator=(const Mat<S_, rows_, cols_> &other) {
    for(int i = 0; i < rows_; i++)
      this->elements[i] = other.elements[i];
    return (*this);
  }
};

template <class S_>
class Mat<S_, 2, 2> : public BasicMat<Mat, S_, 2, 2> {
public:
  Mat() {}

  Mat(S_ s) { this->set(s); }

  Mat(const Mat<S_, 2, 2> &other) {
    (*this) = other;
  }
  
  Mat<S_, 2, 2> &operator=(const Mat<S_, 2, 2> &other) {
    this->elements[0] = other.elements[0];
    this->elements[1] = other.elements[1];
    return (*this);
  }
  
  S_ det() const {
    return (this->elements[0][0] * this->elements[1][1] -
            this->elements[0][1] * this->elements[1][0]);
  }

  Mat<S_, 2, 2> inverse() const {
    Mat<S_, 2, 2> cof_trans;
    const S_ tol = 1E-20;

    // Create matrix of cofactors
    cof_trans[0][0] = this->elements[1][1];
    cof_trans[0][1] = -this->elements[0][1];
    cof_trans[1][0] = -this->elements[1][0];
    cof_trans[1][1] = this->elements[0][0];

    S_ d = det();
    if(fabs(d) < tol) std::cout << "det 2x2 = " << d << std::endl;
    //assert(fabs(d) > tol);

    return cof_trans * (1 / d);
  }
};

template <class S_>
class Mat<S_, 3, 3> : public BasicMat<Mat, S_, 3, 3> {
public:
  Mat() {}

  Mat(S_ s) { this->set(s); }

  Mat(const Mat<S_, 3, 3> &other) {
    (*this) = other;
  }
  
  Mat<S_, 3, 3> &operator=(const Mat<S_, 3, 3> &other) {
    this->elements[0] = other.elements[0];
    this->elements[1] = other.elements[1];
    this->elements[2] = other.elements[2];
    return (*this);
  }
  
  S_ det() const {
    return
      (this->elements[0][0] * this->elements[1][1] * this->elements[2][2] +
       this->elements[0][1] * this->elements[1][2] * this->elements[2][0] +
       this->elements[0][2] * this->elements[1][0] * this->elements[2][1] -
       this->elements[0][2] * this->elements[1][1] * this->elements[2][0] -
       this->elements[0][0] * this->elements[1][2] * this->elements[2][1] -
       this->elements[0][1] * this->elements[1][0] * this->elements[2][2]);
  }

  Mat<S_, 3, 3> inverse() const {
    Mat<S_, 3, 3> cof;
    const S_ tol = 1E-20;

    // Create matrix of cofactors
    for(int i = 0; i < 3; i++)
      for(int j = 0; j < 3; j++) {
        S_ c = this->subMatrix(i, j).det();
        //assert(fabs(c) > tol);
        cof[i][j] = c * ((i + j) % 2 ? -1 : 1);
      }
    
    S_ d = det();
    if(fabs(d) < tol) std::cout << "det 3x3 = " << d << std::endl;
    //assert(fabs(d) > tol);

    return cof.transpose() / d;
  }
};

typedef Mat<float, 2, 2> Mat2f;
typedef Mat<float, 3, 3> Mat3f;

typedef Mat<double, 2, 2> Mat2d;
typedef Mat<double, 3, 3> Mat3d;

typedef Mat<float, 3, 4> Mat3x4f;
typedef Mat<double, 3, 4> Mat3x4d;

#endif
