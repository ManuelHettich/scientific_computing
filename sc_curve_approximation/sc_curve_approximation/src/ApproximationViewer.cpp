//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "ApproximationViewer.h"
#include <imgui.h>
#include "lu.h"
#include "cholesky.h"
#include <Eigen/Dense>
#include <sstream>

//== IMPLEMENTATION ==========================================================

ApproximationViewer::ApproximationViewer(const char* _title, int _width,
                                         int _height)
    : InterpolationViewer(_title, _width, _height)
{
    // set the default values
    fitting_              = INTERPOLATE;
    interpolation_solver_ = LU_EIGEN;
    approximation_solver_ = CHOLESKY_EIGEN;
    poly_degree_          = 0;
}

//-----------------------------------------------------------------------------

void ApproximationViewer::keyboard(int key, int code, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
        return;

    switch (key)
    {
        // switch solvers
        case GLFW_KEY_SPACE:
        {
            if (fitting_ == INTERPOLATE)
            {
                switch (interpolation_solver_)
                {
                    case LU_EIGEN:
                        interpolation_solver_ = LU;
                        break;
                    case LU:
                        interpolation_solver_ = LU_EIGEN;
                        break;
                }
            }
            else
            {
                switch (approximation_solver_)
                {
                    case CHOLESKY_EIGEN:
                        approximation_solver_ = CHOLESKY;
                        break;
                    case CHOLESKY:
                        approximation_solver_ = QR_EIGEN;
                        break;
                    case QR_EIGEN:
                        approximation_solver_ = CHOLESKY_EIGEN;
                        break;
                }
            }

            fit_curve();
            break;
        }

        // simple test case to get a simple decomposition
        case GLFW_KEY_1:
        {
            // get constraints by sampling a sine-curve
            constraints_x_.clear();
            constraints_y_.clear();

            constraints_x_.push_back(-0.5);
            constraints_x_.push_back(0);
            constraints_x_.push_back(0.5);

            constraints_y_.push_back(-1 + 0.5);
            constraints_y_.push_back(0.5);
            constraints_y_.push_back(-1 + 0.5);

            // interpolate constraints
            fit_curve();
            break;
        }

        // test case: interpolate or approximate some sin-curve with different numbers of points
        case GLFW_KEY_2:
        case GLFW_KEY_3:
        case GLFW_KEY_4:
        {
            // key defines the step size dx in x direction
            Scalar dx;
            if (key == GLFW_KEY_2)
                dx = 0.1;
            else if (key == GLFW_KEY_3)
                dx = 0.05;
            else
                dx = 0.025;

            // get constraints by sampling a sine-curve
            constraints_x_.clear();
            constraints_y_.clear();
            for (Scalar x = -1.0; x <= 1.0; x += dx)
            {
                constraints_x_.push_back(x);
                constraints_y_.push_back(x * 0.7 * sin(10.0 * x));
            }

            // interpolate or approximate constraints
            fit_curve();
            break;
        }
        case GLFW_KEY_5:
        case GLFW_KEY_6:
        case GLFW_KEY_7:
        {
            // key defines the step size dx in x direction
            Scalar dx;
            if (key == GLFW_KEY_5)
                dx = 0.1;
            else if (key == GLFW_KEY_6)
                dx = 0.05;
            else
                dx = 0.025;

            // get constraints by sampling a sine-curve
            constraints_x_.clear();
            constraints_y_.clear();
            for (Scalar x = -1.0; x <= 1.0; x += dx)
            {
                constraints_x_.push_back(x);
                constraints_y_.push_back(0.5 * x + 0.25 * sin(7.0 * x));
            }


            // interpolate or approximate constraints
            fit_curve();
            break;
        }

        // c -> clear points and reset curve
        case GLFW_KEY_C:
        {
            coefficients_.clear();
            constraints_x_.clear();
            constraints_y_.clear();
            poly_degree_ = 0;
            break;
        }

        // increase/decrease polynomial degree for approximation
        case GLFW_KEY_F2:
        case GLFW_KEY_F1:
        {
            poly_degree_ =
                (key == GLFW_KEY_F2 ? std::min(poly_degree_+1, (int)constraints_x_.size()-1) : std::max(0, poly_degree_-1));
            std::cout << "Approximate with degree " << poly_degree_
                      << " polynomial\n";
            fit_curve();
            break;
        }

        // switch to approximation mode
        case GLFW_KEY_A:
        {
            fitting_ = APPROXIMATE;
            fit_curve();
            break;
        }

        // switch to interpolation mode
        case GLFW_KEY_I:
        {
            fitting_ = INTERPOLATE;
            fit_curve();
            break;
        }

        // add some noise to the constraints
        case GLFW_KEY_N:
        {
            add_noise(0.01);
            fit_curve();
            break;
        }

        // other keys: let parent handle the key event
        default:
        {
            InterpolationViewer::keyboard(key, code, action, mods);
        }
    }
}

//-----------------------------------------------------------------------------

void ApproximationViewer::processImGUI()
{
    if (ImGui::CollapsingHeader("Fitting", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // interpolation or approximation
        int fitting = (int)fitting_;
        ImGui::RadioButton("Interpolate", &fitting, 0);
        ImGui::RadioButton("Approximate", &fitting, 1);
        if (fitting != fitting_)
        {
            fitting_ = (Fitting_type)fitting;
            fit_curve();
        }

        ImGui::Spacing();
        ImGui::Spacing();
    }

    if(fitting_ == INTERPOLATE)
        InterpolationViewer::processImGUI();
    else
    {
        if (ImGui::CollapsingHeader("Approximation", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // which solver
            int solver = (int)approximation_solver_;
            ImGui::RadioButton("Eigen's Cholesky", &solver, 0);
            ImGui::RadioButton("Our Cholesky",     &solver, 1);
            ImGui::RadioButton("Eigen's QR",       &solver, 2);
            if (solver != approximation_solver_)
            {
                approximation_solver_ = (Solver)solver;
                fit_curve();
            }

            ImGui::Spacing();
            ImGui::Spacing();

            ImGui::PushItemWidth(100);
            poly_degree_ = std::max(0, std::min(poly_degree_, (int)constraints_x_.size()-1));
            int degree = poly_degree_;
            ImGui::SliderInt("Degree", &degree, 0, std::max(0, (int)constraints_x_.size()-1));
            ImGui::PopItemWidth();
            if (degree != poly_degree_)
            {
                poly_degree_ = degree;
                fit_curve();
            }

            ImGui::Spacing();
            ImGui::Spacing();
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::Button("Add Noise")){
        add_noise(0.01);
        fit_curve();
    }
}

//-----------------------------------------------------------------------------

void ApproximationViewer::fit_curve()
{
    if (fitting_ == INTERPOLATE)
        InterpolationViewer::fit_curve();
    else
        approximate();
}

//-----------------------------------------------------------------------------

void ApproximationViewer::add_noise(Scalar _amplitude)
{
/**
 * Add noise to the `y`-component of each constraint.
 * The offset should be between `-_amplitude` and `_amplitude`.
 */

for (auto &y : constraints_y_) {
  y += -_amplitude + static_cast <float> (rand()) /
                         (static_cast <float> (RAND_MAX / (2 * _amplitude)));
}

}

//-----------------------------------------------------------------------------

void ApproximationViewer::approximate()
{
  /**
   * Determine the polynomial coefficients as the
   * least squares solution of the overdetermined problem `A * x = b`.
   * - Setup the overdetermined system `A * x = b`. Make sure to set `m` and `n` correctly.
   *   The polynom's intended degree is stored in `poly_degree_`.
   * - Setup the normal equations.
   * - Depending on the variable `approximation_solver_`, use Eigen's or your own
   *   Cholesky solver to solve the normal equations, or use Eigens's QR
   *   solver to solve the overdetermined system `A * x = b` directly.
   * - Copy the solution in `x` to the `coefficients_` vector.
   *
   * To find out how to use Eigen's Cholesky and QR solvers, simply run a Google search for the terms.
   */

  assert(constraints_x_.size() == constraints_y_.size());

  unsigned int i, j;
  unsigned int m = constraints_x_.size();
  unsigned int n = poly_degree_ + 1;
  MatrixXX A(m, n);
  VectorX b(m), x(n);

  // Fill the system's right hand side 'b'
  for (i = 0; i < m; i++) {
      b(i) = constraints_y_[i];
  }

  // Setup overdetermined system A * x = b
  // Fill matrix 'A'
  for (i = 0; i < m; i++) {
      for (j = 0; j < n ; j++) {
          A(i, j) = pow(constraints_x_[i], j);
      }
  }

  // Setup normal equations
  MatrixXX AtA = A.transpose() * A;
  VectorX Atb = A.transpose() * b;

  // Use correct solver
  switch(approximation_solver_) {

      case CHOLESKY_EIGEN:
      {
          // Eigen Cholesky solver
          x = AtA.ldlt().solve(Atb);

          std::cout << "Error AtA: " << (AtA*x-Atb).norm() << std::endl;

          break;
      }
      case QR_EIGEN:
      {
          // Eigen QR solver
          x = AtA.fullPivHouseholderQr().solve(Atb);

          std::cout << "Error AtA: " << (AtA*x-Atb).norm() << std::endl;

          break;
      }

      case CHOLESKY:
      {
          CholeskySolver solver;
          if (solver.factorize(AtA)) solver.solve(Atb, x);
          break;
      }
  }

  std::cout << "Error A: " << (A*x-b).norm() << std::endl;

  // copy solution to coefficients vector
  coefficients_.resize(n);
  for (i = 0; i < n; ++i) coefficients_[i] = x(i);
}

//=============================================================================
