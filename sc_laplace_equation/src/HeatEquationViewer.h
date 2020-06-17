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

#include <pmp/Window.h>
#include <pmp/MatVec.h>

#include <vector>
#include "types.h"
#include "sparse_matrix.h"

using namespace pmp;


//== CLASS DEFINITION =========================================================


/// A 3D viewer for the heat equation exercise, based on GLFW and ImGUI
class HeatEquationViewer : public Window
{
public:

    /// constructor
    HeatEquationViewer(const char* _title, int _width, int _height);

    /// destructor
    virtual ~HeatEquationViewer();

protected: // GUI functions

    /// this function is called when the scene has to be rendered.
    virtual void display() override;

    /// this function handles keyboard events
    virtual void keyboard(int key, int code, int action, int mods) override;

    /// this function handles mouse button events
    virtual void mouse(int button, int action, int mods) override;

    /// this function handles mouse motion
    virtual void motion(double xpos, double ypos) override;

    /// this function triggers the time_integration() function
    virtual void doProcessing() override;

    /// render/handle GUI
    virtual void processImGUI() override;

    /// pick a point
    void pick_point();

protected:

    /// initialize OpenGL stuff
    void init();

    /// generate a NxN grid
    void generate_grid(int _resolution);

    /// update the normal vectors of the grid vertices.
    /// this function has to be called after modifying the grid values.
    void update_normals();

    /// read-write access grid points
    vec3& grid_point(int i, int j) { return grid_points_[i*grid_resolution_+j]; }

    /// read-write access normal vectors of grid points
    vec3& grid_normal(int i, int j) { return grid_normals_[i*grid_resolution_+j]; }

    /// do one step of time integration using one of the below methods
    void time_integration();

    /// read-only access to x-coordinate of grid point (i,j)
    const float X(int i, int j) { return grid_points_[i*grid_resolution_+j][0]; }

    /// read-only access to y-coordinate of grid point (i,j)
    const float Y(int i, int j) { return grid_points_[i*grid_resolution_+j][1]; }

    /// read-write access to z-coordinate of grid point (i,j)
    float& U(int i, int j) { return grid_points_[i*grid_resolution_+j][2]; }

    /// read-write access to velocities of grid vertex (i,j)
    float& V(int i, int j) { return grid_velocities_[i*grid_resolution_+j]; }

    /// read-write access to acceleration of grid vertex (i,j)
    float& A(int i, int j) { return grid_accelerations_[i*grid_resolution_+j]; }

    /// explicit Euler integration
    void explicit_euler_step();

    /// solve for equilibrium
    void solve_equilibrium();

    /// setup dense system and solve by Eigen's cholesky
    void solve_cholesky();

    /// setup sparse system
    void setup_sparse_system(utils::SparseMatrix& A, utils::Vector& b, utils::Vector x);

    /// solve by our gradient descent solver
    void solve_gd();

    /// solve by our conjugate gradients solver
    void solve_cg();

protected:

    /// grid resolution
    int grid_resolution_;

    /// grid data
    std::vector<vec3>   grid_points_;
    std::vector<vec3>   grid_normals_;
    std::vector<GLuint> grid_indices_;

    /// dynamic data for grid vertices
    std::vector<float>  grid_velocities_, grid_accelerations_;

    /// render setting
    bool render_wireframe_;

    /// color texture
    GLuint textureID_;

    /// which mouse button is active
    int mouse_button_;

    /// switch animation on/off
    bool animate_;

    /// value of time-step
    float time_step_;

    /// time counter
    float integration_time_;

    /// flag to store if left mouse button is being pressed
    bool button_down_;

    /// which solver to use?
    enum Solver
    {
        CHOLESKY_EIGEN=0,
        GRADIENT_DESCENT=1,
        CONJUGATE_GRADIENTS=2
    } solver_;
};

//=============================================================================
