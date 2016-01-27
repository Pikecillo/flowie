/*
 * University of Houston
 * Mario Rincon-Nigro. March 2013.
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "linear_solver.hpp"

//extern "C" {
//#include <taucs.h>
//}

extern "C" void dgesv_(const int *N, const int *nrhs, double *A,
		       const int *lda, int *ipiv, double *b,
		       const int *ldb, int *info);
extern "C" void dgels_(const char *trans, const int *m,
		       const int *n, const int *nrhs, double* a,
		       const int *lda, double *b,
		       const int *ldb, double *work,
		       int *lwork, int *info);
extern "C" void dsyev_(const char *jobz, const char *uplo,
		       const int *n, double *a, const int *lda,
		       double *w, double *work, int *lwork, int *info);

/*extern "C" int taucs_linsolve(taucs_ccs_matrix* A, 
			      void**            F,
			      int               nrhs,
			      void*             X,
			      void*             B,
			      char*             options[],
			      void*             opt_arg[]);*/

/*
 * Computes the solution to a real system of linear equations Ax = B,
 * where A is an N-by-N matrix, and b and x are N-by-NRHS matrices.
 * The LU decomposition with partial pivoting and row interchanges is
 * used to factor A as
 *    A = P * L * U,
 * where P is a permutation matrix, L is unit lower triangular, and U is
 * upper triangular.  The factored form of A is then used to solve the
 * system of equations A * x = b.
 */
void Solver::dgesv(const Matrixd &A, const Matrixd &b, Matrixd &x) {
  int n = A.rows();
  int nrhs = b.cols();
  int lda = n;
  int ldb = n;
  int info;
  int *ipiv;

  // Enforce correct dimensions
  assert(A.rows() == A.cols() && A.rows() == x.rows() &&
	 A.rows() == b.rows() && x.cols() == b.cols());
  
  Matrixd b_copy(b);

  ipiv = (int *)malloc(n * sizeof(int)); 
  
  // Solve the system
  dgesv_(&n, &nrhs, A.getDataPtr(), &lda, ipiv,
	 b_copy.getDataPtr(), &ldb, &info);        
  free(ipiv);
  
  // Make sure lapack did well, otherwise print error
  if(info != 0) {
    std::cout << "error: dgesv info=" << info << std::endl;
    assert(false);
  }

  // Copy solution to x
  for(int i = 0; i < x.rows(); i++)
    for(int j = 0; j < x.cols(); j++) {
      x.set(i, j, b_copy.get(i, j));
    }
}

/*
 * Solves overdetermined or underdetermined real linear
 * systems  involving  an  M-by-N  matrix  A, or its transpose,
 * using a QR or LQ factorization of A. It is assumed A has full rank
 */
void Solver::dgels(const Matrixd &A, const Matrixd &b, Matrixd &x) {
  int m = A.rows();
  int n = A.cols();
  int nrhs = b.cols();
  int lda = m;
  int ldb = (m > n ? m : n);
  double wkopt;
  double *work;
  int lwork;
  int info;
  const char trans[] = "No transpose";
  
  // Enforce correct dimensions
  assert(A.rows() == b.rows() && A.cols() == x.rows() && 
	 x.cols() == b.cols());
  
  Matrixd b_copy(b);
  
  // Query and allocate the optimal workspace
  lwork = -1;
  dgels_(trans, &m, &n, &nrhs,
	 A.getDataPtr(), &lda,
	 b_copy.getDataPtr(), &ldb,
	 &wkopt, &lwork, &info);
  
  lwork = (int)wkopt;
  work = (double *)malloc(lwork * sizeof(double));
  
  // Solve the system
  dgels_(trans, &m, &n, &nrhs, A.getDataPtr(), &lda,
	 b_copy.getDataPtr(), &ldb, work, &lwork, &info);
  free(work);
  
  // Make sure lapack did well
  assert(info == 0);

  // Copy solution to x
  for(int i = 0; i < x.rows(); i++)
    for(int j = 0; j < x.cols(); j++) {
      x.set(i, j, b_copy.get(i, j));
    }
}

void Solver::dsyev(const Matrixd &A,
		   std::vector<Matrixd> &eigenvecs,
		   std::vector<double> &eigenvals) {
  int n = A.cols();
  int lda = n, info, lwork;
  double wkopt;
  double *work;
  double *w;

  // Enforce correct dimensions
  assert(A.rows() == A.cols() && A.cols() > 0);

  /* Local arrays */
  w = (double *)malloc(n * sizeof(double));
  
  // Query and allocate the optimal workspace
  lwork = -1;
  dsyev_("Vectors", "Upper", &n, A.getDataPtr(),
	 &lda, w, &wkopt, &lwork, &info);

  lwork = (int)wkopt;
  work = (double*)malloc(lwork * sizeof(double));
  
  // Solve eigenproblem
  dsyev_("Vectors", "Upper", &n, A.getDataPtr(),
	 &lda, w, work, &lwork, &info);
  
  // Make sure LAPACK did well
  assert(info == 0);

  // Copy eigenvalues
  eigenvals.clear();
  for(int i = 0; i < n; i++)
    eigenvals.push_back(w[i]);

  // Copy eigenvectors
  eigenvecs.clear();
  for(int j = 0; j < n; j++) {
    Matrixd eigvec(n, 1);

    for(int i = 0; i < n; i++)
      eigvec.set(i, 0, A.get(i, j));

    eigenvecs.push_back(eigvec);
  }

  free(work);
}

/*
 * Solves sparse linear system of equations
 */
void Solver::linsolve(const SparseMatrixd &A, const Matrixd &b, Matrixd &x) {
  /*// Create a TAUCS matrix
  taucs_ccs_matrix *A_taucs;

  A_taucs = taucs_dccs_create(A.rows(), A.cols(), A.nonZero());
  A_taucs->flags = TAUCS_DOUBLE;
  
  if(A.isSymmetric())
    A_taucs->flags |= (TAUCS_SYMMETRIC | TAUCS_LOWER);
  else {
    // SparseMatrix is by rows, it is necessary to change it by columns
    SparseMatrixd tmp(A.rows(), A.cols());

  }
  
  // Convert SparseMatrix to the TAUCS CCS
  
  // Solve the linear system
  void *F = NULL;
  const char *options[] = {(const char *)"taucs.factor.LLT=true", NULL };
  void *opt_arg[] = { NULL };

  //taucs_logfile("stdout");

  int i = taucs_linsolve(A_taucs, &F, b.cols(), x.getDataPtr(), b.getDataPtr(),
	   options, opt_arg);
  // Check that TAUCS did well
  if(i != TAUCS_SUCCESS) {
    if(i == TAUCS_ERROR)
      std::cout << "Generic error." << std::endl;
    if(i == TAUCS_ERROR_NOMEM)
      std::cout << "NOMEM error." << std::endl;
    if(i == TAUCS_ERROR_BADARGS)
      std::cout << "BADARGS error." << std::endl;
    if(i == TAUCS_ERROR_MAXDEPTH)
      std::cout << "MAXDEPTH error." << std::endl;
    if(i == TAUCS_ERROR_INDEFINITE)
      std::cout << "NOT POSITIVE DEFINITE error." << std::endl;
  
    assert(false);
    }

  // 5) clean up
  taucs_linsolve(NULL, &F, 0, NULL, NULL, NULL, NULL);
  taucs_dccs_free(A_taucs);*/
}
