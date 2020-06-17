//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "InterpolationViewer.h"
#include <imgui.h>

//== IMPLEMENTATION ==========================================================

InterpolationViewer::
InterpolationViewer(const char* _title, int _width, int _height)
    : pmp::Window(_title, _width, _height)
{
    // start with Eigen's solver
    interpolation_solver_ = LU_EIGEN;

    // OpenGL state
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glDisable( GL_DITHER );
    glDisable( GL_DEPTH_TEST );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}

//-----------------------------------------------------------------------------

void InterpolationViewer::keyboard(int key, int code, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
        return;

    switch (key)
    {
        // space key -> switch between Eigen solver and our solver
        case GLFW_KEY_SPACE:
        {
            if (interpolation_solver_ == LU_EIGEN)
            {
                interpolation_solver_ = LU;
            }
            else
            {
                interpolation_solver_ = LU_EIGEN;
            }
            fit_curve();
            break;
        }


        // c -> clear points and reset curve
        case GLFW_KEY_C:
        {
            coefficients_.clear();
            constraints_x_.clear();
            constraints_y_.clear();
            break;
        }

        // toggle GUI
		case GLFW_KEY_G:
		{
            showImGUI( !showImGUI() );
			break;
		}

        // escape key -> quit
        case GLFW_KEY_Q:
        case GLFW_KEY_ESCAPE:
        {
            exit(0);
            break;
        }
    }
}

//-----------------------------------------------------------------------------

void InterpolationViewer::processImGUI()
{
    if (ImGui::CollapsingHeader("Interpolation", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // which solver
        int solver = (int)interpolation_solver_;
        ImGui::RadioButton("Eigen's LU", &solver, 0);
        ImGui::RadioButton("Our LU",     &solver, 1);
        if (solver != interpolation_solver_)
        {
            interpolation_solver_ = (Solver)solver;
            fit_curve();
        }
    }
}

//-----------------------------------------------------------------------------

void InterpolationViewer::mouse(int button, int action, int mods)
{
    // only on mouse press with left button
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        // query current OpenGL viewing settings
        GLdouble projection[16], modelview[16];
        GLint viewport[4];
        glGetDoublev(GL_PROJECTION_MATRIX, projection);
        glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
        glGetIntegerv(GL_VIEWPORT, viewport);

        // get point under mouse cursor
        double x, y;
        cursorPos(x,y);

        // gluUnproject gives us the point on the near plane
        GLdouble ox, oy, oz;
        gluUnProject(x, m_height-y, 0.0, modelview, projection, viewport, &ox, &oy, &oz);

        if(ox <= 1.0 && ox >= -1.0 && oy < 1.1 && oy > -1.1)
        {
            // add point to interpolation constraints
            constraints_x_.push_back(ox);
            constraints_y_.push_back(oy);

            // re-compute curve
            fit_curve();
        }
    }
}

//-----------------------------------------------------------------------------

void InterpolationViewer::display()
{
    // clear window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // adjust viewport
    int gui_width = 180, width_wo_gui = m_width - gui_width;
    glViewport(gui_width, 0, width_wo_gui, m_height);

	// adjust projection
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	if (width_wo_gui > m_height)
        gluOrtho2D(-1.1 * float(width_wo_gui)/float(m_height),
                    1.1 * float(width_wo_gui)/float(m_height),
                   -1.1,
                    1.1);
    else
        gluOrtho2D(-1.1f,
                    1.1f,
                   -1.1f * float(m_height)/float(width_wo_gui),
                    1.1f * float(m_height)/float(width_wo_gui));

    // draw coordiante axes
    glColor3f(0, 0, 0);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glVertex2f(-1, 0);
    glVertex2f(1, 0);
    glVertex2f(0, -1);
    glVertex2f(0, 1);
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex2f(1.0, 0.01);
    glVertex2f(1.0, -0.01);
    glVertex2f(1.02, 0.0);
    glVertex2f(0.01, 1.0);
    glVertex2f(-0.01, 1.0);
    glVertex2f(0.0, 1.02);
    glEnd();

    // draw interpolation points
    glPointSize(10.0);
    glColor3f(1, 0, 0);
    glBegin(GL_POINTS);
    for (unsigned int i = 0; i < constraints_x_.size(); ++i)
    {
        glVertex2f(constraints_x_[i], constraints_y_[i]);
    }
    glEnd();

    // draw curve (if it has been computed already)
    if (!coefficients_.empty())
    {
        glColor3f(0, 0, 1);
        glLineWidth(3.0);
        glBegin(GL_LINE_STRIP);
        for (Scalar x = -1.0; x <= 1.01; x += 0.01)
        {
            glVertex2f(x, evaluate_curve(x));
        }
        glEnd();
    }
}

//-----------------------------------------------------------------------------

Scalar InterpolationViewer::evaluate_curve(Scalar _x) const
{
    // evalute polynomial as
    // sum_i(coeff[i]*x^i) = coeff[0] + coeff[1]*x + coeff[2]*x^2 + ...

    Scalar value = 0.0;

    for (unsigned int i = 0; i < coefficients_.size(); ++i)
        value += coefficients_[i] * pow(_x, (int)i);

    return value;
}

//-----------------------------------------------------------------------------

void InterpolationViewer::fit_curve()
{
    assert(constraints_x_.size() == constraints_y_.size());

    unsigned int i, j, m = constraints_x_.size();
    MatrixXX A(m, m);
    VectorX b(m), x(m);

    /**
     * Build and solve the system.
     * - Fill matrix `A`. Write `A`'s elements via `A(i,j) = ...`.
     * - Setup the system's right hand side `b`. Write `b`'s elements via `b(i) =
     * ...`.
     * - Use Eigen's LU solver to solve the system `A * x = b`.
     *   See the function `fullPivLu()` of the matrix `A`.
     * - Compute and print the error `|| A * x - b ||`.
     * - Depending on `interpolation_solver_`, use our own LU solver instead.
     */

    // Fill the system's right hand side 'b'
    for (i = 0; i < m; i++) {
      b(i) = constraints_y_[i];
    }

    // Fill matrix 'A'
    for (i = 0; i < m; i++) {
      for (j = 0; j < m ; j++) {
        A(i, j) = pow(constraints_x_[i], j);
      }
    }

    if (interpolation_solver_ == LU_EIGEN) {
      // Use Eigen's LU solver to solve the system 'A * x = b'
      // https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU.html#af563471f6f3283fd10779ef02dd0b748
      x = A.fullPivLu().solve(b);

      /*
      Eigen::FullPivLU<MatrixXX> solver;
      // LU decomp
      solver.compute(A);

      // solve
      x = solver.solve(b);
      */

      std::cout << "Error: " << (A*x-b).norm() << std::endl;

    } else {
      LU_Solver solver_;

      if (solver_.factorize(A)) {
        solver_.solve(b, x);
      }
    }

    // copy solution to coefficients vector
    coefficients_.resize(m);
    for (i = 0; i < m; ++i) coefficients_[i] = x(i);
}

//=============================================================================
