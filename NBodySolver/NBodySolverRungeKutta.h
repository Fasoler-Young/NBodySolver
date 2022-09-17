#pragma once
#include "NBodySolver.h"
#include "NBodeSolverEuler.h"

class NBodySolverRungeKutta :
    public NBodySolverEuler
{
private:
    // Первые нули, для соотвествия номеров
    value_type B2[2] = {0, 2. / 9. };
    value_type B3[3] = {0, 1. / 12., 1. / 4. };
    value_type B4[4] = {0, 69. / 128., -243. / 128., 135. / 64. };
    value_type B5[5] = {0, -17. / 12., 27. / 4., -27. / 5., 16. / 15. };
    value_type B6[6] = {0, 65. / 432., -5. / 16., 13. / 16., 4. / 27., 5 / 144 };

    value_type CH[7] = {0., 47. / 450., 0., 12. / 25., 32. / 225., 1. / 30., 6. / 25. };
    value_type CT[7] = {0., -1. / 150., 0., 3. / 100., -16. / 75., -1. / 20., 6. / 25. };

    //value_type h_new;
    // Здесь k_c 
    std::vector<vector3> k1_c, k2_c, k3_c, k4_c, k5_c, k6_c;
    std::vector<vector3> k1_v, k2_v, k3_v, k4_v, k5_v, k6_v;
    value_type local_err = 1e-18;
    std::vector<vector3> tmp;

    void resize_k(size_t new_size);


public:
    NBodySolverRungeKutta(NBodyData* data);
    void step(value_type* dt);
    void stepRK45(value_type* dt);
    void stepRK4(value_type* dt, vector3* d_coord, vector3* d_v);
    std::string method_name();
    void set_max_local_err(value_type err);

};

