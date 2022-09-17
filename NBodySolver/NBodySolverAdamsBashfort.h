#pragma once
#include "NBodySolver.h"
#include "NBodySolverRungeKutta.h"
class NBodySolverAdamsBashfort :
    public NBodySolver
{
private:
    size_t rank = 5;
    value_type B[5][5] = {
        {1.},
        {3. / 2.,     -1. / 2.},
        {23. / 12.,   -4. / 3.,       5. / 12.},
        {55. / 24.,   -59. / 24.,     37. / 24.,    -3. / 8.},
        {1901. / 720., 1387. / 360.,  109. / 30.,   -637. / 360.,   251. / 720.} };
    NBodySolverRungeKutta* prestep_solver;
    std::vector<std::vector<vector3>> prev_dv, prev_dc;
    std::vector<vector3> corr_c, corr_v;

    void resize(size_t count);
public:
    NBodySolverAdamsBashfort(NBodyData* data);
    void step(value_type* dt);
    std::string method_name();
};

