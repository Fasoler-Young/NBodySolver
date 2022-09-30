#pragma once
#include "NBodeSolverEuler.h"

NBodySolverEuler::NBodySolverEuler(NBodyData* data) : NBodySolver(data)
{
}

void NBodySolverEuler::step(value_type *dt)
{
	vector3* coord = get_data()->get_coord();
	vector3* velosites = get_data()->get_velosites();
	const value_type* mass = get_data()->get_mass();
	size_t count = get_data()->get_count();

	if (!dv.size()) {
		dv.resize(count);
		correction_velosites.resize(count);
		correction_coord.resize(count);
	}

	omp_set_num_threads(7);
	get_data()->calculate_total_force(coord, dv.data());
	#pragma omp parallel for
	for (int id_body1 = 0; id_body1 < count; id_body1++) {
		//const vector3& coord_body1(coord[id_body1]);
		//vector3 total_force;
		//vector3 correction;
		//for (size_t id_body2 = 0; id_body2 < count; id_body2++) {
		//	
		//	if (id_body1 == id_body2) continue;
		//	const vector3& coord_body2(coord[id_body2]);
		//	vector3 force(get_data()->force(coord_body1, coord_body2, mass[id_body1], mass[id_body2]));
		//	//total_force += force;
		//	total_force = get_data()->Kahan_sum(total_force, force, &correction);
		//}
		//dv[id_body1] = total_force / mass[id_body1];
		velosites[id_body1] = get_data()->Kahan_sum(velosites[id_body1], dv[id_body1] * *dt, &correction_velosites[id_body1]);
		
	}
	for (int body = 0; body < count; body++) {
		coord[body] = get_data()->Kahan_sum(coord[body], velosites[body] * *dt, &correction_coord[body]);

	}
	dv.clear();
	get_data()->increase_time(*dt);

}

std::string NBodySolverEuler::method_name()
{
	return "euler";
}
