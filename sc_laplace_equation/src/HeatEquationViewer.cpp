//=============================================================================
//
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//
//=============================================================================

#include "HeatEquationViewer.h"
#include "rainbow.h"
#include "conjugate_gradients.h"
#include "gradient_descent.h"
#include <imgui.h>
#include <pmp/Timer.h>
#include <Eigen/Dense>

using namespace pmp;

//== IMPLEMENTATION ===========================================================

HeatEquationViewer::HeatEquationViewer(const char* _title, int _width, int _height)
    : Window(_title, _width, _height)
{
    // initialize OpenGL stuff
    init();

    // set default values
    render_wireframe_ = true;
    animate_     = false;
    time_step_   = 0.5;
    button_down_ = false;
    integration_time_ = 0.0;
    solver_ = CHOLESKY_EIGEN;

    // initialize grid
    generate_grid(50);
}

//-----------------------------------------------------------------------------

HeatEquationViewer::~HeatEquationViewer()
{
    // clean up texture memory
    if (glIsTexture(textureID_))
        glDeleteTextures( 1, &textureID_);
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::init()
{
    // OpenGL state
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glColor4f(0.0, 0.0, 0.0, 1.0);
    glDisable( GL_DITHER );
    glEnable( GL_DEPTH_TEST );


    // some performance settings
    glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE );
    glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE );


    // surface material
    GLfloat mat_a[] = {0.2, 0.2, 0.2, 1.0};
    GLfloat mat_d[] = {0.8, 0.8, 0.8, 1.0};
    GLfloat mat_s[] = {0.9, 0.9, 0.9, 1.0};
    GLfloat shine[] = {128.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_a);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_d);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_s);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);


    // lighting
    glLoadIdentity();

    GLfloat pos1[] = { 0.1, 0.1, -0.02, 0.0};
    GLfloat pos2[] = {-0.1, 0.1, -0.02, 0.0};
    GLfloat pos3[] = { 0.0, 0.0, 0.1, 0.0};
    GLfloat col1[] = {.05, .05, .6, 1.0};
    GLfloat col2[] = {.6, .05, .05, 1.0};
    GLfloat col3[] = {1.0, 1.0, 1.0, 1.0};

    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_POSITION, pos1);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,  col1);
    glLightfv(GL_LIGHT0,GL_SPECULAR, col3);

    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT1,GL_POSITION, pos2);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,  col2);
    glLightfv(GL_LIGHT1,GL_SPECULAR, col3);

    glEnable(GL_LIGHT2);
    glLightfv(GL_LIGHT2,GL_POSITION, pos3);
    glLightfv(GL_LIGHT2,GL_DIFFUSE,  col3);
    glLightfv(GL_LIGHT2,GL_SPECULAR, col3);

    glEnable(GL_LIGHTING);


    // generate texture
    glGenTextures(1, &textureID_);
    glBindTexture(GL_TEXTURE_2D, textureID_);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, rainbow);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);


    // setup automatic texture coordinates
    glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
    GLfloat plane[4] = { 0.0, 0.0, 5.0, 0.0 };
    glTexGenfv(GL_S, GL_OBJECT_PLANE, plane);


    // setup camera position
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.5, -1.0, 1.0,  // eye
              0.5, 0.5, 0.0,   // center
              0.0, 0.0, 1.0);  // up
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::generate_grid(int _resolution)
{
    // adjust size of grid
    grid_resolution_ = _resolution;
    float grid_spacing = 1.0 / (grid_resolution_-1);


    std::cout
        << "Construct "
        << grid_resolution_ << "x" << grid_resolution_
        << " grid..." << std::flush;


    // allocate memory for points and normal vectors
    grid_points_.resize( grid_resolution_*grid_resolution_ );
    grid_normals_.resize( grid_resolution_*grid_resolution_ );


    // (x,y) values are from unit square [0,1]x[0,1]
    int    i, j;
    float  x, y;
    for (i=0, x=0.0; i<grid_resolution_; ++i, x+=grid_spacing)
    {
        for (j=0, y=0.0; j<grid_resolution_; ++j, y+=grid_spacing)
        {
            grid_point(i,j)  = vec3(x, y, 0.0);
            grid_normal(i,j) = vec3(0.0, 0.0, 1.0);
        }
    }


    // compute indices for OpenGL vertex array rendering
    grid_indices_.clear();
    grid_indices_.reserve((grid_resolution_-1)*(grid_resolution_-1)*4);
    for (i=0; i<grid_resolution_-1; ++i)
    {
        for (j=0; j<grid_resolution_-1; ++j)
        {
            grid_indices_.push_back( (i+0)*grid_resolution_ + (j+0) );
            grid_indices_.push_back( (i+1)*grid_resolution_ + (j+0) );
            grid_indices_.push_back( (i+1)*grid_resolution_ + (j+1) );
            grid_indices_.push_back( (i+0)*grid_resolution_ + (j+1) );
        }
    }


    std::cout << "done\n";


    // reset velocities V and accelerations A
    grid_velocities_.resize( grid_resolution_ * grid_resolution_ );
    grid_accelerations_.resize( grid_resolution_ * grid_resolution_ );
    for (int i=0; i<grid_resolution_; ++i)
        for (int j=0; j<grid_resolution_; ++j)
            V(i,j) = A(i,j) = 0.0;
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::update_normals()
{
    vec3 dx, dy;

    for (int i=0; i<grid_resolution_; ++i)
    {
        for (int j=0; j<grid_resolution_; ++j)
        {
            // compute derivative in x
            if (i-1>=0 && i+1<grid_resolution_)
                dx = grid_point(i+1,j) - grid_point(i-1,j);  // central difference
            else if (i-1>=0)
                dx = grid_point(i,j) - grid_point(i-1,j);    // backward difference
            else
                dx = grid_point(i+1,j) - grid_point(i,j);    // forward difference


            // compute derivative in y
            if (j-1>=0 && j+1<grid_resolution_)
                dy = grid_point(i,j+1) - grid_point(i,j-1);  // central difference
            else if (j-1>=0)
                dy = grid_point(i,j) - grid_point(i,j-1);    // backward difference
            else
                dy = grid_point(i,j+1) - grid_point(i,j);    // forward difference


            // normal is cross product of tangents
            grid_normal(i,j) = normalize(cross(dx, dy));
        }
    }
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::time_integration()
{
    // set up time counters
    static float accumulated_time = 0.0;
    static int ctr = 0;
    // start timing measurement
    ctr++;
    Timer timer;
    timer.start();


    // to one step of time integration
    explicit_euler_step();


    // stop timer, compute elapsed time
    timer.stop();
    accumulated_time += timer.elapsed();
    if(ctr>= 50)
    {
        integration_time_ = accumulated_time/static_cast<float>(ctr);
        ctr = 0;
        accumulated_time = 0.0;
    }


    // update the grid
    update_normals();
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::doProcessing()
{
    if (animate_)
    {
        time_integration();
    }
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::processImGUI()
{
    if (ImGui::CollapsingHeader("Heat Equation", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // animation checkbox
        ImGui::Checkbox("Animate it!", &animate_);

        ImGui::Spacing();
        ImGui::Spacing();

        ImGui::PushItemWidth(100);
        if (ImGui::Button("One Step"))
        {
            time_integration();
        }
        ImGui::SameLine();
        if (ImGui::Button("Equilibrium"))
        {
            solve_equilibrium();
            update_normals();
        }
        ImGui::PopItemWidth();

        ImGui::Spacing();
        ImGui::Spacing();

        // which solver
        int solver = (int)solver_;
        ImGui::RadioButton("Eigen's Cholesky", &solver, 0);
        ImGui::RadioButton("Our Gradient Descent", &solver, 1);
        ImGui::RadioButton("Our Conjugate Gradients", &solver, 2);
        if (solver != solver_)
        {
            solver_ = (Solver)solver;
        }

        ImGui::Spacing();
        ImGui::Spacing();

        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Time Step", &time_step_, 0.01f, 1.2f, "%.2f", 1.5);
        ImGui::PopItemWidth();

        ImGui::Spacing();
        ImGui::Spacing();

        ImGui::PushItemWidth(100);
        int res = grid_resolution_;
        ImGui::SliderInt("Grid Resolution", &res, 10, 100);
        ImGui::PopItemWidth();
        if (res != grid_resolution_)  generate_grid(res);

        ImGui::Spacing();
        ImGui::Spacing();

        ImGui::Checkbox("Render grid", &render_wireframe_);

        ImGui::Spacing();
        ImGui::Spacing();

        if(animate_)
            ImGui::Text("Integration Time: %.3f ms", integration_time_);
    }
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::keyboard(int key, int code, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
        return;

    switch (key)
    {
        // solve laplace equation for equilibrium
        case GLFW_KEY_SPACE:
        {
            solve_equilibrium();
            update_normals();
            break;
        }

        // rotate on left and right keys
        case GLFW_KEY_LEFT:
        {
            glMatrixMode(GL_MODELVIEW);
            glTranslatef(0.5, 0.5, 0.0);
            glRotatef(-5.0, 0.0, 0.0, 1.0);
            glTranslatef(-0.5, -0.5, 0.0);
            break;
        }
        case GLFW_KEY_RIGHT:
        {
            glMatrixMode(GL_MODELVIEW);
            glTranslatef(0.5, 0.5, 0.0);
            glRotatef(5.0, 0.0, 0.0, 1.0);
            glTranslatef(-0.5, -0.5, 0.0);
            break;
        }

        // generate some test cases
        case GLFW_KEY_1:
        {
            for (int i = 0; i < grid_resolution_; ++i)
                for (int j = 0; j < grid_resolution_; ++j)
                    U(i, j) = 0.20 * pow(sin(3.0 * M_PI * X(i, j)) * sin(3.0 * M_PI * Y(i, j)), 2.0);

            update_normals();
            break;
        }

        case GLFW_KEY_2:
        {
            for (int i = 0; i < grid_resolution_; ++i)
                for (int j = 0; j < grid_resolution_; ++j)
                    U(i, j) =
                        0.1 * X(i, j) + 0.1 * Y(i, j) +
                        0.2 * pow(sin(4.0 * M_PI * X(i, j)) * sin(4.0 * M_PI * Y(i, j)), 3.0);

            update_normals();
            break;
        }

        case GLFW_KEY_3:
        {
            for (int i = 0; i < grid_resolution_; ++i)
                for (int j = 0; j < grid_resolution_; ++j)
                    U(i, j) =
                        0.20 *
                        pow(cos(1.0 * M_PI * X(i, j)) * sin(3.0 * M_PI * Y(i, j)), 2.0);

            update_normals();
            break;
        }

        // toggle automatic animation on/off
        case GLFW_KEY_A:
        {
            animate_ = !animate_;
            break;
        }

        // do just one step of time integration
        case GLFW_KEY_S:
        {
            time_integration();
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

void HeatEquationViewer::mouse(int button, int action, int mods)
{
    // only on mouse press with left button
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        button_down_ = true;
        mouse_button_ = button;
        pick_point();
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
    {
        button_down_ = false;
    }
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::motion(double _x, double _y)
{
    if (button_down_)
    {
        pick_point();
    }
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::pick_point()
{
    // query current OpenGL viewing settings
    GLdouble  modelview[16], projection[16];
    GLint     viewport[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);

    // get point under mouse cursor
    double x, y;
    cursorPos(x,y);

    // reverse projection to get two 3D points = ray
    int    y_new = viewport[3] - y; // in OpenGL y=0 is at the 'bottom'
    double p0[3], p1[3];
    gluUnProject(x, y_new, 0.0,
            modelview,
            projection,
            viewport,
            &p0[0], &p0[1], &p0[2]);
    gluUnProject(x, y_new, 1.0,
            modelview,
            projection,
            viewport,
            &p1[0], &p1[1], &p1[2]);


    // intersect ray with z=0 plane
    float l = -p0[2] / (p1[2]-p0[2]);
    vec3 p(p0[0] + l*(p1[0]-p0[0]),
            p0[1] + l*(p1[1]-p0[1]),
            p0[2] + l*(p1[2]-p0[2]));


    // lift point and some neighbors
    for (int i = 0; i < grid_resolution_; ++i)
    {
        for (int j = 0; j < grid_resolution_; ++j)
        {
            // q is grid point (i,j) with z-coordinate set to zero.
            vec3 q = grid_point(i, j);
            q[2]    = 0.0;

            // if q is close enough to point _p on z=0 plane, update q's z-value
            if (norm(q - p) < 0.05) U(i, j) = 0.2;
        }
    }

    // recompute grid normals
    update_normals();
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::display()
{
    // clear screen
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // adjust viewport
    glViewport(0, 0, m_width, m_height);

    // adjust projection
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective(45.0, (GLfloat)m_width/(GLfloat)m_height, 0.1, 10.0);
    glMatrixMode( GL_MODELVIEW );

    // setup OpenGL vertex arrays
    glVertexPointer(3, GL_FLOAT, 0, (const float*)(grid_points_[0].data()));
    glEnableClientState(GL_VERTEX_ARRAY);
    glNormalPointer(GL_FLOAT, 0, (const float*)(grid_normals_[0].data()));
    glEnableClientState(GL_NORMAL_ARRAY);

    // draw filled faces
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthRange(0.01, 1.0);
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_1D);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_TEXTURE_GEN_S);
    glDrawElements(GL_QUADS, grid_indices_.size(), GL_UNSIGNED_INT, &grid_indices_[0]);
    glDisable(GL_TEXTURE_GEN_S);
    glDisable(GL_TEXTURE_1D);
    glDepthRange(0.0, 1.0);

    // draw grid edges
    if (render_wireframe_)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glColor3f(0,0,0);
        glDisable(GL_LIGHTING);
        glDrawElements(GL_QUADS, grid_indices_.size(), GL_UNSIGNED_INT, &grid_indices_[0]);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    // reset OpenGL state
    glDisable(GL_VERTEX_ARRAY);
    glDisable(GL_NORMAL_ARRAY);
    glDisable(GL_LIGHTING);
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::explicit_euler_step()
{
    /**
     * \todo Update the height values `U(i, j)` at every grid point `(i, j)` by
     * an explicit integration (Euler) of the heat equation.
     *
     * The grid `U` and `V` that contain the height values and the velocities,
     * respectively, already exist. Just take care of the following two steps:
     * 1. Compute the time-derivative of `U(i, j)` and store the value in
     *    the velocity field `V(i, j)` at every point `(i, j)`.
     * 2. Update `U(i, j)` in the direction `V(i, j)` for every point `(i, j)`.
     *
     * Hints:
     * - `time_step_` is the time step, required for updating `U`.
     * - `grid_resolution_` is the grid resolution, i.e. both `i` and `j` run
     *   from `0` to `grid_resolution - 1`. Points at the boundary should stay fixed.
     * - Consider the grid spacing to be `2`.
     */
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::solve_equilibrium()
{
    // solve equilibrium
    switch (solver_)
    {
        case CHOLESKY_EIGEN:
        {
            solve_cholesky();
            break;
        }

        case GRADIENT_DESCENT:
        {
            solve_gd();
            break;
        }

        case CONJUGATE_GRADIENTS:
        {
            solve_cg();
            break;
        }
    }
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::solve_cholesky()
{
    Timer timer;
    timer.start();

    /**
     * \todo Update the height values `U(i, j)` at every grid point `(i, j)` by
     * solving the Laplace equation. This will compute the equilibrium of
     * the heat flow.
     *
     * 1. Setup the linear system `A * x = b` for the Laplace equation.
     * 2. Since the system is symmetric positive definite, use a Cholesky
     * solver, provided by Eigen, to solve it.
     *
     * Hints:
     * - Again, boundary points stay fixed, so the matrix `A` is
     *   `(grid_resolution_ - 2) * (grid_resolution_ - 2)` by
     *   `(grid_resolution_ - 2) * (grid_resolution_ - 2)`.
     * - Take care when mapping a grid index `(i, j)` to an index `k` in the vector `x`
     *   and the matrix `A`. Use the mapping `k = idx(i, j)` discussed in the lecture.
     * - You can ignore the grid spacing `h`, or set it to `1`, since the factor
     *   `1 / h^2` will cancel out.
     * - Note that height values of boundary grid points have to be moved to the
     *   right hand side `b`. Otherwise, the system cannot be solved.
     */

    timer.stop();
    std::cerr << "Cholesky Solve: " << timer << std::endl;
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::setup_sparse_system(utils::SparseMatrix& A, utils::Vector& b, utils::Vector x)
{
    const int r = grid_resolution_;
    int i, j, k;

#define IDX(i, j) ((i)-1) * (r - 2) + ((j)-1)

    // setup A & b
    for (i = 1; i < r - 1; ++i)
    {
        for (j = 1; j < r - 1; ++j)
        {
            A.begin_row();

            k = IDX(i, j);

            // initialize x
            x(k) = U(i, j);

            // initialize right hand side
            b(k) = 0.0;

            // diagonal
            A.add_value(k, 4.0);

            // left neighbor
            if (i > 1)
                A.add_value(IDX(i - 1, j), -1.0);
            else
                b(k) += U(i - 1, j);

            // right neighbor
            if (i < r - 2)
                A.add_value(IDX(i + 1, j), -1.0);
            else
                b(k) += U(i + 1, j);

            // bottom neighbor
            if (j > 1)
                A.add_value(IDX(i, j - 1), -1.0);
            else
                b(k) += U(i, j - 1);

            // top neighbor
            if (j < r - 2)
                A.add_value(IDX(i, j + 1), -1.0);
            else
                b(k) += U(i, j + 1);

            A.end_row();
        }
    }

#undef IDX
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::solve_gd()
{
    Timer timer;
    timer.start();

    const int r = grid_resolution_;
    const int m = (r - 2) * (r - 2);
    int i, j;
    utils::SparseMatrix A(m, m);
    utils::Vector x(m), b(m);

    // setup matrix A, right-hand side b, and initial guess for x
    setup_sparse_system(A, b, x);

    // solve using gradient descent, up to error 1e-10, with up to 1000 iterations
    gradient_descent(A, b, x, 1e-10, 1000);

    // copy solution x to grid U(i,j)
#define IDX(i, j) ((i)-1) * (r - 2) + ((j)-1)
    for (i = 1; i < r - 1; ++i)
        for (j = 1; j < r - 1; ++j) U(i, j) = x(IDX(i, j));
#undef IDX

    // output timing & error
    timer.stop();
    std::cout << "Gradient Descent:\n";
    std::cout << "  time:  " << timer << std::endl;
    Scalar err = (A * x - b).norm() / b.norm();
    std::cout << "  error: " << err << " --> "
              << (err < 1e-10 ? "ok\n" : "TOO HIGH\n");
}

//-----------------------------------------------------------------------------

void HeatEquationViewer::solve_cg()
{
    Timer timer;
    timer.start();

    const int r = grid_resolution_;
    const int m = (r - 2) * (r - 2);
    int i, j;
    utils::SparseMatrix A(m, m);
    utils::Vector x(m), b(m);

    // setup matrix A, right-hand side b, and initial guess for x
    setup_sparse_system(A, b, x);

    // solve using conjugate gradients, up to error 1e-10, with up to 1000
    // iterations
    conjugate_gradients(A, b, x, 1e-10, 1000);

    // copy solution x to grid U(i,j)
#define IDX(i, j) ((i)-1) * (r - 2) + ((j)-1)
    for (i = 1; i < r - 1; ++i)
        for (j = 1; j < r - 1; ++j) U(i, j) = x(IDX(i, j));
#undef IDX

    // output timing & error
    timer.stop();
    std::cout << "Conjugate Gradients:\n";
    std::cout << "  time:  " << timer << std::endl;
    Scalar err = (A * x - b).norm() / b.norm();
    std::cout << "  error: " << err << " --> "
              << (err < 1e-10 ? "ok\n" : "TOO HIGH\n");
}

//=============================================================================
