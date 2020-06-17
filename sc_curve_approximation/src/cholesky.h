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

#include <iostream>
#include <Eigen/Dense>
#include "lu.h"

//== CLASS DEFINITION =========================================================

/// Our Cholesky solver
class CholeskySolver : public LU_Solver
{
public:

    /// empty constructor
    CholeskySolver() {}

    /// factorize matrix A=L*L^T
    virtual bool factorize(const MatrixXX& _A) override;

    /// solve A*x=b
    virtual void solve(const VectorX& _b, VectorX& _x) override;
};

//=============================================================================
