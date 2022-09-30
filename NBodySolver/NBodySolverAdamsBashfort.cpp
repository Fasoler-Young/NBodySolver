#include "NBodySolverAdamsBashfort.h"

void NBodySolverAdamsBashfort::resize(size_t count)
{
	prev_dc.resize(rank);
	prev_dv.resize(rank);
	for (size_t i = 0; i < rank; i++) {
		prev_dv[i].resize(count);
		prev_dc[i].resize(count);
	}
	corr_c.resize(count);
	corr_v.resize(count);
}

NBodySolverAdamsBashfort::NBodySolverAdamsBashfort(NBodyData* data) : NBodySolver(data)
{
	prestep_solver = new NBodySolverRungeKutta(data);
}

void NBodySolverAdamsBashfort::step(value_type* dt)
{
	vector3* coord = get_data()->get_coord();
	vector3* velosites = get_data()->get_velosites();
	const value_type* mass = get_data()->get_mass();
	size_t count = get_data()->get_count();

	if (prev_dc.size() != rank)
		resize(count);
	size_t step = get_data()->get_step();
	if (step < rank) {
		prestep_solver->stepRK4(dt, prev_dc[step].data(), prev_dv[step].data());
	}
	else {
		// y_n+k = y_n+k-1 + sum_1_k(B_k*dt*f(y_n+k))
		// ћы будем хранить в предыдущих значени€х dt*f(y_n+k), потому что rk4 так делает и потому что так проще
		//std::vector<vector3> dv(count), dcoord(count);
		prev_dc.erase(prev_dc.begin());
		prev_dv.erase(prev_dv.begin());
		prev_dc.push_back(std::vector<vector3> (count));
		prev_dv.push_back(std::vector<vector3>(count));
		get_data()->calculate_total_force(coord, prev_dv[prev_dv.size() - 1].data());

#pragma omp parallel for
		for (int body = 0; body < count; body++) {
			prev_dc[prev_dc.size() - 1][body] = velosites[body] * *dt;
			vector3 dcoord;
			for (size_t i = 0; i < rank; i++) {
				dcoord += prev_dc[rank - i - 1][body] * B[rank - 1][i];
			}
			coord[body] = get_data()->Kahan_sum(coord[body], dcoord, &corr_c[body]);
		}
#pragma omp parallel for
		for (int body = 0; body < count; body++) {
			prev_dv[prev_dv.size() - 1][body] *= *dt;
			vector3 dv;
			for (size_t i = 0; i < rank; i++) {
					dv += prev_dv[rank - i - 1][body] * B[rank - 1][i];
			}
			velosites[body] = get_data()->Kahan_sum(velosites[body], dv, &corr_v[body]);
		}


		get_data()->increase_time(*dt);
	}
}

std::string NBodySolverAdamsBashfort::method_name()
{
	return "adams";
}
