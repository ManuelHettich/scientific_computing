//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================
#pragma once
//=============================================================================

#include "sc_math.h"
#include "sparse_matrix.h"
#include <fstream>

//== CLASS DEFINITION =========================================================

///  Gradient Descent Solver
void gradient_descent(const utils::SparseMatrix& A, const utils::Vector& b, utils::Vector& x,
                      Scalar _error, int _iters)
{
    assert(A.rows() == A.cols());

    /**
     * \todo Implement the gradient descent solver as presented in the lecture.
     *
     * Besides the function parameters `x` and `b`, you will have to store only
     * the residual `r`, and the product `Ar = A * r`.
     *
     * The stopping criterion should be the relative residual error, i.e.,
     * `norm(A * x - b) / norm(b) < _error`. Note that `r = b - A * x` is computed anyway,
     * as well as `r * r`, so you don't have to compute `(A * x - b)` again for the
     * error check.
     */
}

//=============================================================================
