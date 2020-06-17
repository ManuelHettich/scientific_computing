//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "cholesky.h"

#include <Eigen/Dense>

//== CLASS DEFINITION =========================================================

bool CholeskySolver::factorize(const MatrixXX& A)
{
  assert(A.rows() == A.cols());

  int i, j, k;
  const int m = A.rows();

  // initialize L
  L.resize(m, m);
  L.setZero();

  /**
   * Compute the Cholesky factorization, i.e., compute the matrix `L`.
   * The matrix `A` is given as a parameter, matrix `L` is a member variable
   * (inherited from LU_Solver). The matrices can be accessed by `A(i, j)` and
   * `L(i, j)`. The function should return `true` on success and `false` on
   * failure.
   */

  // Initialise L with A
  L = A;

  // set upper triangle of L to zero
  for (i = 0; i < m; i++) {
    for (j = i + 1; j < m; j++) {
      L (i, j) = 0.0;
    }
  }

  for (k = 0; k < m; k++)
  {
    const Scalar diag = L(k, k);

    /*
     * Check for numerical stability:
     * fabs(diag) should be > eps for 1/diag
     * diag should be > 0 for sqrt(diag)
     */
    if (fabs(diag) < 5 * std::numeric_limits<Scalar>::min()) {
      std::cerr << "CholeskySolver: Factorization failed.\n";
        return false;
    }

    // rank-one update of L(k:m, k:m) in lower triangle
    for (i = k + 1; i < m; i++)
      for (j = i; j < m; j++) L(j, i) -= L(j, k) * (L(i, k) / diag);

    // divide column k by sqrt of diagonal
    for (i = k; i < m; i++) L(i, k) /= sqrt(diag);
  }

  // check error of factorization
  std::cout << "  error(A=L*L^T) : " << (A - L * L.transpose()).norm()
            << std::endl;

  return true;
}

//-----------------------------------------------------------------------------

void CholeskySolver::solve(const VectorX& _b, VectorX& _x)
{
  // use LU's solve function with U = L^T
  // Note: If memory was limited, this would not be recommended.
  U = L.transpose();
  LU_Solver::solve(_b, _x);
}

//=============================================================================
