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
#include "lu.h"
#include <vector>


//== CLASS DEFINITION =========================================================

/// Simple 2D viewer for the first exercise, based on GLFW and ImGUI
class InterpolationViewer : public pmp::Window
{
public:

    /// constructor
    InterpolationViewer(const char* _title, int _width, int _height);

protected:

    /// render/handle GUI
    virtual void processImGUI() override;

    /// this function handles keyboard events
    virtual void keyboard(int key, int code, int action, int mods) override;

    /// this function handles mouse button events
    virtual void mouse(int button, int action, int mods) override;

    /// this function is called when the scene has to be rendered.
    virtual void display() override;

    /// fit a polynomial curve to constraints
    virtual void fit_curve();

    /// evaluate the polynomial specified by coefficients_
    Scalar evaluate_curve(Scalar _x) const;

protected:

    /// which solver to use?
    enum Solver { LU_EIGEN=0, LU=1 } interpolation_solver_;

    /// x-coordinates of constraints
    std::vector<Scalar> constraints_x_;

    /// y-coordinates of constraints
    std::vector<Scalar> constraints_y_;

    /// coefficients of the polynomial
    std::vector<Scalar> coefficients_;
};

//=============================================================================
