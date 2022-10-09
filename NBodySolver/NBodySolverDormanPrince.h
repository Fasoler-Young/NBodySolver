#pragma once
#include"NBodySolver.h"

class NBodySolverDormanPrince :
    public NBodySolver
{
public:
    std::vector<std::vector<value_type>> B;

    // 4 и 5 пор€док соотвественно
    std::vector<value_type> B1, B2;
    //value_type B1[7] = { 35. / 384.,        0.,     500. / 1113.,    125. / 192.,  -2187. / 6784.,    11. / 84.,      0. };
    //value_type B2[7] = { 5179. / 57600.,    0.,     7571. / 16695.,  393. / 640.,  -92097. / 339200., 187. / 2100.,   1. / 40. };

    // ѕервый индекс - тело, второй - ранг
    std::vector<std::vector<vector3>> k_c, k_v;
    // Ќа задаче двух тел при local_err = 1e-22 отстутствует рост ошибки (при сравнивании координат)
    value_type local_err = 1e-18;
    size_t rank = 7;
    value_type rounding_error;
    bool adaptive_step = true;

    void resize_k(size_t count);
    value_type calculate_err(NBodyData* data, value_type dt);
    bool need_restep(value_type err);
    value_type calculate_new_dt(value_type err, value_type dt);
    value_type max_4(value_type a, value_type b, value_type c, value_type d);


    NBodySolverDormanPrince(NBodyData* data);
    void step(value_type* dt);

    void set_max_local_err(value_type err);
    std::string method_name();
};


class RK4:public NBodySolverDormanPrince {
private:

public:
    RK4(NBodyData* data):NBodySolverDormanPrince(data){
        rank = 4;

        B = {
               {0.},
               {0.5},
               {0., 0.5},
               {0., 0., 1.} };

        B2 = { 1. / 6., 1. / 3., 1. / 3., 1. / 6. };
        adaptive_step = false;
    }
    std::string method_name() { return "rk4"; }
};
