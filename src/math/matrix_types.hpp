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

#ifndef __MATRIX_TYPES_HPP__
#define __MATRIX_TYPES_HPP__

#include <math.h>

#include <iostream>

#include <math/vector_types.hpp>

namespace flowie
{

template <template <typename, int, int> class M, typename S,
	  int Rows, int Cols>
class BasicMat {
protected:
    void set(S s) {
	for(int i = 0; i < Rows; i++)
	    elements[i] = Vector<S, Cols>(s);
    }
    
public:
    const Vector<S, Cols> operator[](int i) const {
	assert(i >= 0 && i <= Rows);
	return elements[i];
    }
    
    Vector<S, Cols> &operator[](int i) {
	assert(i >= 0 && i <= Rows);
	return elements[i];
    }
    
    const Vector<S, Cols> row(int i) const {
	assert(i >= 0 && i <= Rows);
	return elements[i];
    }
    
    Vector<S, Cols> &row(int i) {
	assert(i >= 0 && i <= Rows);
	return elements[i];
    }
    
    Vector<S, Rows> col(int j) const {
	assert(j >= 0 && j <= Cols);
	Vector<S, Rows> c;
	for(int i = 0; i < Rows; i++)
	    c[i] = elements[i][j];
	return c;
    }
    
    M<S, Rows, Cols> &operator+=(const M<S, Rows, Cols> &other) {
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		elements[i][j] += other.elements[i][j];
	
	return *((M<S, Rows, Cols> *)this);
    }
    
    M<S, Rows, Cols> &operator-=(const M<S, Rows, Cols> &other) {
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		elements[i][j] -= other.elements[i][j];
	
	return *((M<S, Rows, Cols> *)this);
    }
    
    M<S, Rows, Cols> &operator*=(const S &scalar) {
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		elements[i][j] *= scalar;
	
	return *((M<S, Rows, Cols> *)this);
    }
    
    template <int Rowso, int Colso>
    M<S, Rows, Colso> &operator*=(const M<S, Rowso, Colso> &other) {
	M<S, Rows, Colso> result(*((M<S, Rows, Colso> *)this));
	
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Colso; j++) {
		S dot = 0.0;
		
		for(int k = 0; k < Cols; k++)
		    dot += (result.elements[i][k] * other.elements[k][j]);
		
		elements[i][j] = dot;
	    }
	
	return *((M<S, Rows, Colso> *)this); 
    }
    
    M<S, Rows, Cols> &operator/=(const S &scalar) {
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		elements[i][j] /= scalar;
	
	return *((M<S, Rows, Cols> *)this);
    }
    
    M<S, Rows, Cols> operator+(const M<S, Rows, Cols> &other) const {
	M<S, Rows, Cols> result(*((M<S, Rows, Cols> *)this));
	result += other;
	return result;
    }
    
    M<S, Rows, Cols> operator-() const {
	M<S, Rows, Cols> result(*((M<S, Rows, Cols> *)this));
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		result[i][j] = -elements[i][j];
	
	return result;
    }
    
    M<S, Rows, Cols> operator-(const M<S, Rows, Cols> &other) const {
	M<S, Rows, Cols> result(*((M<S, Rows, Cols> *)this));
	result -= other;
	return result;
    }
    
    M<S, Rows, Cols> operator*(const S &scalar) const {
	M<S, Rows, Cols> result(*((M<S, Rows, Cols> *)this));
	result *= scalar;
	return result;
    }
    
    Vector<S, Rows> operator*(const Vector<S, Cols> &vector) const {
	Vector<S, Rows> result;
	for(int i = 0; i < Rows; i++)
	    result[i] = elements[i].dot(vector);
	
	return result;
    }
    
    template <int Rowso, int Colso>
    M<S, Rows, Colso> operator*
    (const M<S, Rowso, Colso> &other) const {
	M<S, Rows, Colso> result(*((M<S, Rows, Colso> *)this));
	result *= other;
	return result;
    }
    
    M<S, Rows, Cols> operator/(const S &scalar) const {
	M<S, Rows, Cols> result(*((M<S, Rows, Cols> *)this));
	result /= scalar;
	return result;
    }
    
    M<S, Cols, Rows> transpose() const {
	M<S, Cols, Rows> trans(*((M<S, Cols, Rows> *)this));
	
	for(int i = 0; i < Cols; i++)
	    for(int j = 0; j < Rows; j++)
		trans[i][j] = this->elements[j][i];
	
	return trans;
    }
    
    M<S, Rows - 1, Cols - 1> subMatrix(int i, int j) const {
	assert(i >= 0 && i < Rows && j >= 0 && j < Cols);
	
	M<S, Rows - 1, Cols - 1> m;
	for(int x = 0; x < Rows - 1; x++)
	    for(int y = 0; y < Cols - 1; y++)
		m[x][y] =
		    this->elements[(x < i ? x : x + 1)][(y < j ? y : y + 1)];
	return m;
    }
    
    static M<S, Rows, Cols> identity() {
	M<S, Rows, Cols> id;
	
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		id[i][j] = (i == j ? (S)1 : (S)0);
	
	return id;
    }
    
    template <int Dim>
    static M<S, Dim, Dim> outerProduct(const Vector<S, Dim> &a,
					 const Vector<S, Dim> &b) {
	M<S, Dim, Dim> outer;
	
	for(int i = 0; i < Rows; i++)
	    for(int j = 0; j < Cols; j++)
		outer[i][j] = a[i] * b[j];
	
	return outer;
    }
    
    template <int Dim>
    static M<S, Dim, Dim> projection(const Vector<S, Dim> &v) {
	return (M<S, Dim, Dim>::identity() -
		M<S, Dim, Dim>::outerProduct(v, v));
    }
    
protected:
    Vector<S, Cols> elements[Rows]; 
};

template <typename S, int Rows, int Cols>
class Mat : public BasicMat<Mat, S, Rows, Cols> {
public:
    Mat() { assert(Rows > 0 && Cols > 0); }
    
    Mat(S s) { this->set(s); }
    
    Mat(const Mat<S, Rows, Cols> &other) {
	(*this) = other;
    }
    
    Mat<S, Rows, Cols> &operator=(const Mat<S, Rows, Cols> &other) {
	for(int i = 0; i < Rows; i++)
	    this->elements[i] = other.elements[i];
	return (*this);
    }
};

template <class S>
class Mat<S, 2, 2> : public BasicMat<Mat, S, 2, 2> {
public:
    Mat() {}
    
    Mat(S s) { this->set(s); }
    
    Mat(const Mat<S, 2, 2> &other) {
	(*this) = other;
    }
    
    Mat<S, 2, 2> &operator=(const Mat<S, 2, 2> &other) {
	this->elements[0] = other.elements[0];
	this->elements[1] = other.elements[1];
	return (*this);
    }
    
    S det() const {
	return (this->elements[0][0] * this->elements[1][1] -
		this->elements[0][1] * this->elements[1][0]);
    }
    
    Mat<S, 2, 2> inverse() const {
	Mat<S, 2, 2> cof_trans;
	const S tol = 1E-20;
	
	cof_trans[0][0] = this->elements[1][1];
	cof_trans[0][1] = -this->elements[0][1];
	cof_trans[1][0] = -this->elements[1][0];
	cof_trans[1][1] = this->elements[0][0];
	
	S d = det();
	assert(fabs(d) > tol);
	
	return cof_trans * (1 / d);
    }
};

template <class S>
class Mat<S, 3, 3> : public BasicMat<Mat, S, 3, 3> {
public:
    Mat() {}
    
    Mat(S s) { this->set(s); }
    
    Mat(const Mat<S, 3, 3> &other) {
	(*this) = other;
    }
    
    Mat<S, 3, 3> &operator=(const Mat<S, 3, 3> &other) {
	this->elements[0] = other.elements[0];
	this->elements[1] = other.elements[1];
	this->elements[2] = other.elements[2];
	return (*this);
    }
    
    S det() const {
	return
	    (this->elements[0][0] *
	     this->elements[1][1] *
	     this->elements[2][2] +
	     this->elements[0][1] *
	     this->elements[1][2] *
	     this->elements[2][0] +
	     this->elements[0][2] *
	     this->elements[1][0] *
	     this->elements[2][1] -
	     this->elements[0][2] *
	     this->elements[1][1] *
	     this->elements[2][0] -
	     this->elements[0][0] *
	     this->elements[1][2] *
	     this->elements[2][1] -
	     this->elements[0][1] *
	     this->elements[1][0] *
	     this->elements[2][2]);
    }
    
    Mat<S, 3, 3> inverse() const {
	Mat<S, 3, 3> cof;
	const S tol = 1E-20;
	
	for(int i = 0; i < 3; i++)
	    for(int j = 0; j < 3; j++) {
		S c = this->subMatrix(i, j).det();
		assert(fabs(c) > tol);
		cof[i][j] = c * ((i + j) % 2 ? -1 : 1);
	    }
	
	S d = det();
	assert(fabs(d) > tol);
	
	return cof.transpose() / d;
    }
};

typedef Mat<float, 2, 2> Mat2f;
typedef Mat<float, 3, 3> Mat3f;

typedef Mat<double, 2, 2> Mat2d;
typedef Mat<double, 3, 3> Mat3d;

typedef Mat<float, 3, 4> Mat3x4f;
typedef Mat<double, 3, 4> Mat3x4d;

}

#endif
