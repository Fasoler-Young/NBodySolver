#include "NBodySolverDormanPrince.h"
#include <iostream>

NBodySolverDormanPrince::NBodySolverDormanPrince(NBodyData* data) :NBodySolver(data)
{
	value_type tmp = 1.;
	while (1 + tmp > 1) {
		rounding_error = tmp;
		tmp /= 2;
	}
}

value_type NBodySolverDormanPrince::calculate_new_dt(value_type err, value_type dt) {
	// h_new = h * min (fac_max, max(fac_min,	fac * (local_err / err)^1/5))
	//		   dt		fac_max		1/fac_max	|			tmp				|
	// fac_max ~ [1.5:5];
	value_type fac_max = 3.; 
	value_type tmp = 0.9 * pow((local_err / err), 1. / 5.);
	value_type max = 1 / fac_max > tmp ? 1 / fac_max : tmp;
	value_type min = fac_max < max ? fac_max : max;
	return dt * min;
}


bool NBodySolverDormanPrince::need_restep(value_type err) {
	return false;
}


value_type NBodySolverDormanPrince::calculate_err(size_t count) {
	//return fabs(get_data()->last_total_energy() - get_data()->total_energy());

	// eps = sqrt(1/n*sum((dt((x_end_i_4 - x_end_i_5)/max(10^-5, |x_end_i_4|, |x_start_i_4|, 2*round_err/local_err))^2)


	value_type err = -666.;
	for (int body = 0; body < count; body++) {
		value_type tmp = 0.;
		for (size_t i = 0; i < rank; i++) {
			tmp += (B2[i] - B1[i]) * k_c[body][i].length();
		}
		if (err < tmp)
			err = tmp;
	}
	return fabs(err);
}


void NBodySolverDormanPrince::resize_k(size_t count) {
	k_c.resize(count);
	k_v.resize(count);
	for (size_t i = 0; i < count; i++) {
		k_c[i].resize(rank);
		k_v[i].resize(rank);
	}
}


void NBodySolverDormanPrince::step(value_type* dt)
{
	vector3* coord = get_data()->get_coord();
	vector3* velosites = get_data()->get_velosites();
	const value_type* mass = get_data()->get_mass();
	size_t count = get_data()->get_count();
	if (k_c.size() == 0)
		resize_k(count);

	value_type k_c_err = 0.;
	value_type k_v_err = 0.;
	value_type old_dt = *dt;

	for (size_t i = 0; i < rank; i++) {
		#pragma omp parallel for
		for (int body = 0; body < count; body++) {
			vector3 d_coord;
			for (size_t j = 0; j < i; j++) {
				d_coord += k_c[body][j] * B[i][j];
			}
			vector3 total_force (get_data()->calculate_total_force(d_coord, body));
			k_v[body][i] = total_force / mass[body] * *dt;
			k_c[body][i] = velosites[body]* *dt;
			for (size_t k = 0; k < i; k++) {
				k_c[body][i] += k_v[body][k] * B[i][k] * *dt;
			}
		}
	}

	if (adaptive_step) {
		value_type err = calculate_err(count);
		*dt = calculate_new_dt(err, *dt);
		if (need_restep(err)) {
			std::cout << "Restep: old_dt: " << old_dt << "\tnew_dt: " << *dt << "\terr: " << err << std::endl;
			step(dt);
			return;
		}
	}
	//#pragma omp parallel for
	for (int body = 0; body < count; body++) {
		for (size_t i = 0; i < rank; i++) {
			velosites[body] += k_v[body][i] * B1[i];
			coord[body] += k_c[body][i] * B1[i];
		}
	}

	get_data()->increase_time(old_dt);

}

void NBodySolverDormanPrince::set_max_local_err(value_type err)
{
	local_err = err;
}

std::string NBodySolverDormanPrince::method_name()
{
	return "rkdp";
}
