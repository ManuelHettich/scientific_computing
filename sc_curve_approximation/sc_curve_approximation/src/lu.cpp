//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "lu.h"

//== CLASS DEFINITION =========================================================

bool LU_Solver::factorize(const MatrixXX& A)
{
  assert(A.rows() == A.cols());

  int i, j, k;
  const int m = A.rows();

  // initialize L, U
  L.resize(m, m);
  L = MatrixXX::Identity(m, m);
  U.resize(m, m);
  U = A;

  // main loop --> subtraction row
  for (k = 0; k < m; ++k)
  {
    // if the next diagonal element is too small, the matrix is singular
    if (fabs(U(k, k)) < 5 * std::numeric_limits<Scalar>::min())
    {
      std::cerr << "LU_Solver: Factorization failed." << std::endl;
      return false;
    }

    // Adjust the submatrix, i.e. write to `L` and `U`

    // start with row k+1
    for (i = k + 1; i < m; i++)
    {
      // 1) Platziere Faktoren von L
      L(i, k) = U(i, k) / U(k, k);
      // 2) Subtrahiere Vielfache der Reihe k von den unteren Reihen
      for (j = k; j < m; j++)
      {
        U(i, j) -= L(i, k) * U(k, j);
      }
    }

    // std::cout << "A: \n" << A << std::endl;
    // std::cout << "L: \n" << L << std::endl;
    // std::cout << "U: \n" << U << std::endl;
  }

  // check error of LU factorization
  std::cout << "  error(A = L*U) : " << (A - L * U).norm() << std::endl;
  return true;
}

//-----------------------------------------------------------------------------

void LU_Solver::solve(const VectorX& _b, VectorX& _x)
{
  /**
   * Solve the system `A * _x = b`, using the computed factorization of
   * the matrix `A = L * U`. The right hand side is `_b`, the result is to be
   * written to `_x`.
   * - Solve `L * y = b`, then check the error `norm(L * y - b)`.
   * - Solve `U * x = y`, then check the error `norm(U * x - y)`.
   */

  int i, j;
  const int m = L.rows();
  VectorX _y(m);
  double sum = 0.0;

  // 1) Solve `L * y = b`
  for (i = 0; i < m; i++)
  {

    sum = 0.0;
    for (j = 0; j < i; j++)
    {
      sum += L(i, j) * _y(j);
    }
    _y(i) = (_b(i) - sum) / L(i, i);
  }

  // std::cout << "  y: \n" << _y << std::endl;
  std::cout << "  error(L * y = b) : " << (L * _y - _b).norm() << std::endl;

  // 2) Solve `U * x = y`
  for (i = m-1; i >= 0; i--)
  {
    sum = 0.0;
    for (j = i + 1; j < m; j++)
    {
      sum += U(i, j) * _x(j);
    }
    _x(i) = (_y(i) - sum) / U(i, i);
  }

  // std::cout << "  x: \n" << _x << std::endl;
  std::cout << "  error(U * x = y) : " << (U * _x - _y).norm() << std::endl;
}

//=============================================================================
