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

#include "InterpolationViewer.h"

//== CLASS DEFINITION =========================================================

/// Simple 2D viewer for the second exercise, based on GLFW and ImGUI
class ApproximationViewer : public InterpolationViewer
{
public:

    /// constructor
    ApproximationViewer(const char* _title, int _width, int _height);

protected:

    /// render/handle GUI
    virtual void processImGUI() override;

    /// this function handles keyboard events
    virtual void keyboard(int key, int code, int action, int mods) override;

    /// add noise to the y-coordinate of constraints
    void add_noise(Scalar _amplitude);

    /// overload this method: interpolate OR approximate constraints
    virtual void fit_curve() override;

    /// compute polynomial that approximates the constraints
    void approximate();

protected:

    /// interpolate or approximate
    enum Fitting_type
    {
        INTERPOLATE=0,
        APPROXIMATE=1
    } fitting_;

    /// which solver to use?
    enum Solver
    {
        CHOLESKY_EIGEN=0,
        CHOLESKY=1,
        QR_EIGEN=2
    } approximation_solver_;

    int poly_degree_;
};

//=============================================================================
