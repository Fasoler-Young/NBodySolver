#include "NBodySolverDormanPrince.h"
#include <iostream>

NBodySolverDormanPrince::NBodySolverDormanPrince(NBodyData* data) :NBodySolver(data)
{
	B = {
		{ 0.},
		{ 1. / 5.},
		{ 3. / 40.,          9. / 40. },
		{ 44. / 45.,         -56. / 15.,      32. / 9 },
		{ 19372. / 6561.,    -25360. / 2187,  64448. / 6561.,    -212 / 729. },
		{ 9017. / 3168.,     -355. / 33.,     46732. / 5247.,    49. / 176.,     -5103. / 18656. },
		{ 35. / 384.,        0.,              500. / 1113.,      125. / 192.,    -2187. / 6784.,     11. / 84. } };
	B2 = { 5179. / 57600.,    0.,     7571. / 16695.,  393. / 640.,  -92097. / 339200., 187. / 2100.,   1. / 40. };
	B1 = { 35. / 384., 0., 500. / 1113., 125. / 192., -2187. / 6784., 11. / 84., 0. };
	value_type tmp = 1.;
	while (1 + tmp > 1) {
		rounding_error = tmp;
		tmp /= 2;
	}
}

value_type NBodySolverDormanPrince::calculate_new_dt(value_type err, value_type dt) {
	// h_new = h/max(0.1, min (5, ((err/local_err)^1/5)/0.9))
	
	value_type tmp = pow((err / local_err), 1. / 5.) / 0.9;
	value_type min = 5 < tmp ? 5 : tmp;
	value_type max = 0.5 > min ? 0.5 : min;
	return dt / max;


	//// h_new = h * min (fac_max, max(fac_min,	fac * (local_err / err)^1/5))
	////		   dt		fac_max		1/fac_max	|			tmp				|
	//// fac_max ~ [1.5:5];

	//value_type fac_max = 3.; 
	//value_type tmp = 0.9 * pow((local_err / err), 1. / 5.);
	//value_type max = 1 / fac_max > tmp ? 1 / fac_max : tmp;
	//value_type min = fac_max < max ? fac_max : max;
	//return dt * min;
}

value_type NBodySolverDormanPrince::max_4(value_type a, value_type b, value_type c, value_type d)
{
	value_type max = a > b ? a : b;
	value_type max2 = c > d ? c : d;
	return max > max2 ? max : max2;
}


bool NBodySolverDormanPrince::need_restep(value_type err) {
	return false;
}


value_type NBodySolverDormanPrince::calculate_err(NBodyData* data, value_type dt) {
	//return fabs(get_data()->last_total_energy() - get_data()->total_energy());

	size_t count = data->get_count();
	vector3* coord = data->get_coord();
	value_type err = 0.;

	for (int body = 0; body < count; body++) {
		vector3 d_coord_4, d_coord_5;
		for (size_t i = 0; i < rank; i++) {
			d_coord_4 += k_c[i][body] * B1[i];
			d_coord_5 += k_c[i][body] * B2[i];
		}

		value_type tmp = (d_coord_4 - d_coord_5).length();
		if (tmp > err)
			err = tmp;
	}
	return err;


	// eps = sqrt(1/n*sum((dt((x_end_i_4 - x_end_i_5)/max(10^-5, |x_end_i_4|,	|x_start_i_4|, 2*round_err/local_err))^2)
	//							d_coord_4   d_coord_5			fabs(d_coord_4)	 coord[body]			
	//size_t count = data->get_count();
	//vector3* coord = data->get_coord();
	//value_type err = 0.;

	//for (int body = 0; body < count; body++) {
	//	vector3 d_coord_4, d_coord_5;
	//	for (size_t i = 0; i < rank; i++) {
	//		d_coord_4 += k_c[i][body] * B1[i];
	//		d_coord_5 += k_c[i][body] * B2[i];
	//	}
	//	vector3 max = vector3(
	//		max_4(1e-5, fabs(d_coord_4.x), fabs(coord[body].x), 2 * rounding_error / local_err),
	//		max_4(1e-5, fabs(d_coord_4.y), fabs(coord[body].y), 2 * rounding_error / local_err),
	//		max_4(1e-5, fabs(d_coord_4.z), fabs(coord[body].z), 2 * rounding_error / local_err)
	//	);
	//	value_type tmp = sqrt((((d_coord_4 - d_coord_5) * dt) / max).norm() / 3.);
	//	if (tmp > err)
	//		err = tmp;
	//}
	//return err;


	//value_type err = -666.;
	//for (int body = 0; body < count; body++) {
	//	value_type tmp = 0.;
	//	for (size_t i = 0; i < rank; i++) {
	//		tmp += (B2[i] - B1[i]) * k_c[body][i].length();
	//	}
	//	if (err < tmp)
	//		err = tmp;
	//}
	//return fabs(err);
}


void NBodySolverDormanPrince::resize_k(size_t count) {
	k_c.resize(rank);
	k_v.resize(rank);
	for (size_t i = 0; i < rank; i++) {
		k_c[i].resize(count);
		k_v[i].resize(count);
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

	std::vector<vector3> tmp_coord(count);


	for (size_t i = 0; i < rank; i++) {
		if(i == 0)
			get_data()->calculate_total_force(coord, k_v[i].data());
		else {
			#pragma omp parallel for
			for (int body = 0; body < count; body++) {
				vector3 d_coord;
				for (size_t j = 0; j < i; j++) {
					d_coord += k_c[j][body] * B[i][j];
				}
				tmp_coord[body] = coord[body] + d_coord;
			}
			get_data()->calculate_total_force(tmp_coord.data(), k_v[i].data());
		}
		#pragma omp parallel for
		for (int body = 0; body < count; body++) {

			//vector3 total_force ();
			//k_v[body][i] = total_force / mass[body] * *dt;
			k_v[i][body] *= *dt;
			k_c[i][body] = velosites[body]* *dt;
			for (size_t k = 0; k < i; k++) {
				k_c[i][body] += k_v[k][body] * B[i][k] * *dt;
			}

		}
	}

	if (adaptive_step) {
		value_type err = calculate_err( get_data(), *dt);
		*dt = calculate_new_dt(err, *dt);
		//std::cout << err << "\t new_dt: " << *dt << std::endl;

		if (need_restep(err)) {
			std::cout << "Restep: old_dt: " << old_dt << "\tnew_dt: " << *dt << "\terr: " << err << std::endl;
			step(dt);
			return;
		}
	}
	#pragma omp parallel for
	for (int body = 0; body < count; body++) {
		vector3 cor_c, cor_v;
		for (size_t i = 0; i < rank; i++) {
			velosites[body] = get_data()->Kahan_sum(velosites[body], k_v[i][body] * B2[i], &cor_v);
			coord[body] = get_data()->Kahan_sum(coord[body], k_c[i][body] * B2[i], &cor_c);
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
