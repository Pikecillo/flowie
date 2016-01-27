/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#ifndef __PCA_HPP__
#define __PCA_HPP__

#include <assert.h>

#include <vector>

#include "linear_solver.hpp"

template <class T_>
T_ mean(std::vector<T_> &data) {
  T_ m = (T_)0;

  for(unsigned int i = 0; i < data.size(); i++)
    m += data[i];

  return m / (data.size());
}

template <class T_>
T_ covariance(std::vector<T_> &data1, std::vector<T_> &data2) {
  assert(data1.size() == data2.size());

  T_ mean1 = mean(data1), mean2 = mean(data2);
  T_ cov = (T_)0;
  int n = data1.size();

  for(int i = 0 ; i < n; i++)
    cov += ((data1[i] - mean1) * (data2[i] - mean2) / (n - 1));

  return cov;
}

template <class T_>
Matrix<T_> covariance_matrix(std::vector<std::vector<T_> > data) {
  Matrix<T_> cov_mat (data.size(), data.size());
  std::vector<T_> avg;

  for(unsigned int i = 0; i < data.size(); i++)
    avg.push_back(mean(data[i]));

  for(unsigned int i = 0; i < data.size(); i++)
    for(unsigned int j = 0; j < data.size(); j++) {
      if(i > j) continue;

      assert(data[i].size() == data[j].size());

      T_ cov = (T_)0;
      int n = data[i].size();

      for(unsigned int k = 0 ; k < data[i].size(); k++)
	cov += ((data[i][k] - avg[i]) * (data[j][k] - avg[j]) / (n - 1));

      cov_mat.at(i, j) = cov;
      cov_mat.at(j, i) = cov;
    }

  return cov_mat;
}

template <class T_>
void pca(std::vector<std::vector<T_> > data) {
  Matrix<T_> cov_mat(data.size(), data.size());
  std::vector<Matrixd> eigenvectors;
  std::vector<double> eigenvalues;

  cov_mat = covariance_matrix(data);

  Solver::dsyev(cov_mat, eigenvectors, eigenvalues);
}

#endif
