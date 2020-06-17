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

int main(int argc, char **argv)
{
    HeatEquationViewer viewer("Heat Equation", 800, 800);
    return viewer.run();
}
