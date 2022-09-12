#pragma once

#include "NBodySolver.h"

class NBodySolverEuler :
    public NBodySolver
{
    std::vector<vector3> dv;
    std::vector<vector3> correction_coord;
    std::vector<vector3> correction_velosites;
public:
    NBodySolverEuler(NBodyData* data);
    void step(value_type *dt);
};

