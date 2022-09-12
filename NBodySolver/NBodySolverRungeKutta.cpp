#pragma once
#include "NBodySolverRungeKutta.h"


#include<iostream>

void NBodySolverRungeKutta::resize_k(size_t new_size)
{
	k1_c.resize(new_size);
	k2_c.resize(new_size);
	k3_c.resize(new_size);
	k4_c.resize(new_size);
	k5_c.resize(new_size);
	k6_c.resize(new_size);
	k1_v.resize(new_size);
	k2_v.resize(new_size);
	k3_v.resize(new_size);
	k4_v.resize(new_size);
	k5_v.resize(new_size);
	k6_v.resize(new_size);

	tmp.resize(new_size);
}

NBodySolverRungeKutta::NBodySolverRungeKutta(NBodyData* data) : NBodySolverEuler(data) {

}

void NBodySolverRungeKutta::step(value_type* dt)
{
	stepRK45(dt);
}

void NBodySolverRungeKutta::stepRK4(value_type* dt)
{
	vector3* coord = get_data()->get_coord();
	vector3* velosites = get_data()->get_velosites();
	const value_type* mass = get_data()->get_mass();
	size_t count = get_data()->get_count();

	if (k1_c.size() != count)
		resize_k(count);

	#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force( get_data()->calculate_total_force(vector3(), body_id));
		k1_v[body_id] = total_force / mass[body_id] * *dt;
		k1_c[body_id] = velosites[body_id] * *dt;
		//tmp[body_id] = coord[body_id] + k1_c[body_id] / 2.;
	}
	#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k1_c[body_id] / 2., body_id));
		k2_v[body_id] = total_force / mass[body_id] * *dt;
		k2_c[body_id] = (velosites[body_id] + k1_v[body_id] / 2.) * *dt;
		//tmp[body_id] = coord[body_id] + k2_c[body_id] / 2.;
	}
	#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k2_c[body_id] / 2., body_id));
		k3_v[body_id] = total_force / mass[body_id] * *dt;
		k3_c[body_id] = (velosites[body_id] + k2_v[body_id] / 2.) * *dt;
		//tmp[body_id] = coord[body_id] + k3_c[body_id] / 2.;
	}
	#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k3_c[body_id] / 2., body_id));
		k4_v[body_id] = total_force / mass[body_id] * *dt;
		k4_c[body_id] = (velosites[body_id] + k3_v[body_id] / 2.) * *dt;

		velosites[body_id] += (k1_v[body_id] + k2_v[body_id] * 2 + k3_v[body_id] * 2 + k4_v[body_id]) / 6.;
		coord[body_id] += (k1_c[body_id] + k2_c[body_id] * 2 + k3_c[body_id] * 2 + k4_c[body_id]) / 6.;
	}
	get_data()->increase_time(*dt);
}




// При eps = 1e-15 рост ошибки отсутствует для задачи двух тел
// Метод показывает достаточную точность за нормальное время. 
// При повышении порядка eps точность растет заметно.
// Показал на этой же задаче точность как у rk4 при eps = 1e-17
void NBodySolverRungeKutta::stepRK45(value_type* dt)
{
	
	vector3* coord = get_data()->get_coord();
	vector3* velosites = get_data()->get_velosites();
	const value_type* mass = get_data()->get_mass();
	size_t count = get_data()->get_count();

	if (k1_c.size() != count)
		resize_k(count);

#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(vector3(), body_id));
		k1_v[body_id] = total_force / mass[body_id] * *dt;
		k1_c[body_id] = velosites[body_id] * *dt;
		//tmp[body_id] = coord[body_id] + k1_c[body_id] / 2.;
	}
#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k1_c[body_id] * B2[1], body_id));
		k2_v[body_id] = total_force / mass[body_id] * *dt;
		k2_c[body_id] = (velosites[body_id] + k1_v[body_id] * B2[1]) * *dt;
		//tmp[body_id] = coord[body_id] + k2_c[body_id] / 2.;
	}
#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k1_c[body_id] * B3[1] + k2_c[body_id] * B3[2], body_id));
		k3_v[body_id] = total_force / mass[body_id] * *dt;
		k3_c[body_id] = (velosites[body_id] + k1_v[body_id] * B3[1] + k2_v[body_id] * B3[2]) * *dt;
		//tmp[body_id] = coord[body_id] + k3_c[body_id] / 2.;
	}
#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k1_c[body_id] * B4[1] + k2_c[body_id] * B4[2] + k3_c[body_id] * B4[3], body_id));
		k4_v[body_id] = total_force / mass[body_id] * *dt;
		k4_c[body_id] = (velosites[body_id] + k1_v[body_id] * B4[1] + k2_v[body_id] * B4[2] + k3_v[body_id] * B4[3]) * *dt;
	}
#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k1_c[body_id] * B5[1] + k2_c[body_id] * B5[2] + k3_c[body_id] * B5[3] + k4_c[body_id] * B5[4], body_id));
		k5_v[body_id] = total_force / mass[body_id] * *dt;
		k5_c[body_id] = (velosites[body_id] + k1_v[body_id] * B5[1] + k2_v[body_id] * B5[2] + k3_v[body_id] * B5[3] + k4_v[body_id] * B5[4]) * *dt;
		//tmp[body_id] = coord[body_id] + k3_c[body_id] / 2.;
	}
	value_type min_h_new = 10;
	value_type max_TE = 0.;
	value_type TE = 0.;
	value_type eps = 1e-15;
#pragma omp parallel for
	for (int body_id = 0; body_id < count; body_id++) {
		vector3 total_force(get_data()->calculate_total_force(k1_c[body_id] * B6[1] + k2_c[body_id] * B6[2] + k3_c[body_id] * B6[3] + k4_c[body_id] * B6[4] + k5_c[body_id] * B6[5], body_id));
		k6_v[body_id] = total_force / mass[body_id] * *dt;
		k6_c[body_id] = (velosites[body_id] + k1_v[body_id] * B6[1] + k2_v[body_id] * B6[2] + k3_v[body_id] * B6[3] + k4_v[body_id] * B6[4] + k5_v[body_id] * B6[5]) * *dt;
		//tmp[body_id] = coord[body_id] + k3_c[body_id] / 2.;
		
		//coord[body_id] += k1_c[body_id] * CH[1] + k2_c[body_id] * CH[2] + k3_c[body_id] * CH[3] + k4_c[body_id] * CH[4] + k5_c[body_id] * CH[5] + k6_c[body_id] * CH[6];
		//velosites[body_id] += k1_v[body_id] * CH[1] + k2_v[body_id] * CH[2] + k3_v[body_id] * CH[3] + k4_v[body_id] * CH[4] + k5_v[body_id] * CH[5] + k6_v[body_id] * CH[6];

		TE = (k1_c[body_id] * CT[1] + k2_c[body_id] * CT[2] + k3_c[body_id] * CT[3] + k4_c[body_id] * CT[4] + k5_c[body_id] * CT[5] + k6_c[body_id] * CT[6]).length();
		value_type temp = 0.9 * *dt * pow((eps / TE), 1. / 5.);

		if (max_TE < TE)
			max_TE = TE;
		if (temp < min_h_new)
			min_h_new = temp;
		//if (TE > 1e-5) {
		//	dt = h_new;
		//	// repeat step with h_new
		//}
	}
	if (max_TE <= eps) {
		#pragma omp parallel for
		for (int body_id = 0; body_id < count; body_id++) {
			coord[body_id] += k1_c[body_id] * CH[1] + k2_c[body_id] * CH[2] + k3_c[body_id] * CH[3] + k4_c[body_id] * CH[4] + k5_c[body_id] * CH[5] + k6_c[body_id] * CH[6];
			velosites[body_id] += k1_v[body_id] * CH[1] + k2_v[body_id] * CH[2] + k3_v[body_id] * CH[3] + k4_v[body_id] * CH[4] + k5_v[body_id] * CH[5] + k6_v[body_id] * CH[6];
		}
		get_data()->increase_time(*dt);
		*dt = min_h_new;
	}
	else {
		std::cout << "Restep: " << *dt << '\t' << min_h_new << '\t' << max_TE << std::endl;
		*dt = min_h_new;
		stepRK45(dt);
	}


	//for (size_t body1 = 0; body1 < get_data()->get_count(); body1++) {
	//	// Здесь мы полагаем, что скорость на протяжении шага постоянна, тогда
	//	// k_v = dt*(F(x)/m)=dt*dv(x), k_c = dt*(dt*dv(x)) = dt*k_v
	//	
	//	//vector3 total_force = get_data()->calculate_total_force(B2[1], body1);
	//	vector3 k1_v(get_data()->calculate_total_force(vector3(), body1) * dt);
	//	vector3 k1_c(k1_v * dt);
	//	// Это может не работать, т.к. смещение по иксу выбирается на основании предыдущего k,
	//	// а они различаются на порядки
	//	// Кажется снова полагая, что скорость постоянна, ориентируемся на k_c
	//	vector3 k2_v(get_data()->calculate_total_force(k1_c * B2[1], body1) * dt);
	//	vector3 k2_c(k2_v * dt);
	//	vector3 k3_v(get_data()->calculate_total_force(k1_c * B3[1] + k2_c * B3[2], body1) * dt);
	//	vector3 k3_c(k3_v * dt);
	//	vector3 k4_v(get_data()->calculate_total_force(k1_c * B4[1] + k2_c * B4[2] + k3_c * B4[3], body1) * dt);
	//	vector3 k4_c(k4_v * dt);
	//	vector3 k5_v(get_data()->calculate_total_force(k1_c * B5[1] + k2_c * B5[2] + k3_c * B5[3] + k4_c * B5[4], body1) * dt);
	//	vector3 k5_c(k5_v * dt);
	//	vector3 k6_v(get_data()->calculate_total_force(k1_c * B6[1] + k2_c * B6[2] + k3_c * B6[3] + k4_c * B6[4] + k5_c * B6[5], body1) * dt);
	//	vector3 k6_c(k6_v * dt);

	//	coord[body1] += k1_c * CH[1] + k2_c * CH[2] + k3_c * CH[3] + k4_c * CH[4] + k5_c * CH[5] + k6_c * CH[6];
	//	velosites[body1] += k1_v * CH[1] + k2_v * CH[2] + k3_v * CH[3] + k4_v * CH[4] + k5_v * CH[5] + k6_v * CH[6];

	//}
	//get_data()->increase_time(dt);

	
}
