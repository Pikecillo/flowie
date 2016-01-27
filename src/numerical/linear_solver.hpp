/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __LINSOLVER_HPP__
#define __LINSOLVER_HPP__

#include <assert.h>

#include <map>
#include <vector>

/*
 * Linear Systems Solvers
 */

/*
 * Matrix is stored in column-major order
 */
template <class T>
class Matrix {
public:
  Matrix(int r, int c) : m_rows(r), m_cols(c) {
    assert(r * c > 0);

    m_data = new T[m_rows * m_cols];
  }
  
  Matrix(int r, int c, T *data) : m_rows(r), m_cols(c) {
    assert(r * c > 0);

    m_data = new T[m_rows * m_cols];
    
    for(int i = 0; i < m_rows * m_cols; i++)
      m_data[i] = data[i];
  }
  
  Matrix(const Matrix<T> &other) {
    m_cols = other.m_cols;
    m_rows = other.m_rows;
    m_data = new T[m_cols * m_rows];
    
    for(int i = 0; i < m_rows * m_cols; i++)
      m_data[i] = other.m_data[i];
  }
  
  ~Matrix() {
    delete [] m_data;
  }
  
  int cols() const { return m_cols; }
  int rows() const { return m_rows; }
  
  T &at(int i, int j) {
    assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
    return m_data[j * m_rows + i];
  }
  
  void set(int i, int j, T value) {
    assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
    m_data[j * m_rows + i] = value;
  }
  
  T get(int i, int j) const {
    assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
    return m_data[j * m_rows + i];
  }
  
  T *getDataPtr() const {
    return m_data;
  }
  
private:
  int m_rows, m_cols;
  T *m_data;
};

template <class T>
class SparseMatrix {
public:
  SparseMatrix(int r, int c, bool s=false)
    : symmetric(s), ncols(c), nnZero(0), data(r) {}

  T get(int i, int j) const {
    assert(i >= 0 && i < data.size() && j >= 0 && j < ncols);

    return data[i][j];
  }

  void set(int i, int j, T value) const {
    assert(i >= 0 && i < data.size() && j >= 0 && j < ncols);

    if(data[i].count(j) == 0) nonZero++;
    data[i][j] = value;
  }

  std::map<int, T> row(int i) const {
    assert(i >= 0 && i < data.size());
    return data[i];
  }

  int cols() const { return ncols; }
  int rows() const { return data.size(); }
  int nonZero() const { return nnZero; }
  bool isSymmetric() const { return symmetric; }

private:
  bool symmetric;
  int ncols;
  int nnZero;
  std::vector<std::map<int, T> > data;
};

typedef Matrix<double> Matrixd;
typedef Matrix<float> Matrixf;

typedef SparseMatrix<double> SparseMatrixd;

class Solver {
public:
  static void dgesv(const Matrixd &A, const Matrixd &b, Matrixd &x);
  static void dgels(const Matrixd &A, const Matrixd &b, Matrixd &x);
  static void dsyev(const Matrixd &A,
		    std::vector<Matrixd> &eigenvecs,
		    std::vector<double> &eigenvals);
  static void linsolve(const SparseMatrixd &A, const Matrixd &b, Matrixd &x);
};

#endif
