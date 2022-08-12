#pragma once

#include "NBodySolver.h"

class NBodeSolverEuler :
    private NBodySolver
{
    std::vector<vector3> dv;
    std::vector<vector3> correction_coord;
    std::vector<vector3> correction_velosites;
public:
    NBodeSolverEuler(NBodyData* data);
    void step(value_type dt);
};

