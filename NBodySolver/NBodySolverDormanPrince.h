#pragma once
#include"NBodySolver.h"

class NBodySolverDormanPrince :
    public NBodySolver
{
private:
    value_type B[7][7] = { 
        { 0.},
        { 1. / 5.},
        { 3. / 40.,          9. / 40. },
        { 44. / 45.,         -56. / 15.,      32. / 9 },
        { 19372. / 6561.,    -25360. / 2187,  64448. / 6561.,    -212 / 729. },
        { 9017. / 3168.,     -355. / 33.,     46732. / 5247.,    49. / 176.,     -5103. / 18656. },
        { 35. / 384.,        0.,              500. / 1113.,      125. / 192.,    -2187. / 6784.,     11. / 84. } };

    // 4 и 5 пор€док соотвественно
    value_type B1[7] = { 35. / 384.,        0.,     500. / 1113.,    125. / 192.,  -2187. / 6784.,    11. / 84.,      0. };
    value_type B2[7] = { 5179. / 57600.,    0.,     7571. / 16695.,  393. / 640.,  -92097. / 339200., 187. / 2100.,   1. / 40. };

    // ѕервый индекс - тело, второй - ранг
    std::vector<std::vector<vector3>> k_c, k_v;
    // Ќа задаче двух тел при local_err = 1e-22 отстутствует рост ошибки (при сравнивании координат)
    value_type local_err = 1e-18;
    size_t rank = 7;
    value_type rounding_error;
    bool adaptive_step = true;

    void resize_k(size_t count);
    value_type calculate_err(size_t count);
    bool need_restep(value_type err);
    value_type calculate_new_dt(value_type err, value_type dt);

public:
    NBodySolverDormanPrince(NBodyData* data);
    void step(value_type* dt);

    void set_max_local_err(value_type err);
    std::string method_name();
};

