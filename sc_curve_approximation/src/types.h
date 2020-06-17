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

#include "Eigen/Dense"

//=============================================================================

/// we use double as Scalar type
//typedef float Scalar;
typedef double Scalar;

typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXX;     ///< A n by m matrix.
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;                   ///< A nd column vector.

//=============================================================================
