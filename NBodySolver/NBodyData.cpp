#include "NBodyData.h"

NBodyData::NBodyData()
{
	time = 0;
	step = 0;
	count = 0;

	current_total_potential_energy = 0.;
	prev_tottal_potential_energy = 0.;
	current_total_kinetic_energy = 0.;
	prev_total_kinetic_energy = 0.;
	current_total_impulse = vector3();
	prev_total_impulse = vector3();


}

vector3 NBodyData::force(const vector3& coord1, const vector3& coord2, const value_type mass1, const value_type mass2)
{
	vector3 dr(coord1 - coord2);
	value_type r2 = dr.norm();
	if (r2 < NBODY_MIN_RADIUS) {
		r2 = NBODY_MIN_RADIUS;
	}
	value_type scalar = (-G * mass1 * mass2 / (r2 * sqrt(r2)));
	return dr*scalar;
}

void NBodyData::increase_time(value_type dt)
{
	time += dt;
	step++;
}

// ≈сли радиус меньше минимального взаимодействие все равно есть, так что, 
// возможно, следует и потенциальную энергию высчитывать
value_type NBodyData::potential_energy(size_t body1, size_t body2) const
{
	value_type dr = (coord[body1] - coord[body2]).length();
	if ( dr < NBODY_MIN_RADIUS) 
		return dr = NBODY_MIN_RADIUS;
	return -G * mass[body1] * mass[body2] / dr;
}

value_type NBodyData::get_time() const
{
	return time;
}

size_t NBodyData::get_step() const
{
	return step;
}

vector3* NBodyData::get_coord()
{
	return coord.data();
}

const vector3* NBodyData::get_coord() const
{
	return coord.data();
}

vector3* NBodyData::get_velosites()
{
	return velosites.data();
}

const vector3* NBodyData::get_velosites() const
{
	return velosites.data();
}

const value_type* NBodyData::get_mass() const
{
	return mass.data();
}

void NBodyData::load_galaxy(const char* path)
{
	auto loc = std::locale(std::locale(), new MyCharCType());
	filestream.imbue(loc);

	filestream.ignore(std::numeric_limits<std::streamsize>::max(),
		filestream.widen('\n'));
	filestream.open(path, std::ios::in);
	value_type data1;
	std::vector<value_type> param;
	while (filestream >> data1)
	{
		param.push_back(data1);
		if (param.size() == 7) {
			add_body(param);
			param.clear();
		}
	}
	filestream.close();
}

void NBodyData::add_body(std::vector<value_type> param)
{
	coord.push_back(vector3(param[0], param[1], param[2]));
	velosites.push_back(vector3(param[3], param[4], param[5]));
	mass.push_back(param[6]);
	count++;
}

void NBodyData::add_body(vector3 coord, vector3 velosites, value_type mass)
{
	this->coord.push_back(coord);
	this->velosites.push_back(velosites);
	this->mass.push_back(mass);
	this->count++;
}

void NBodyData::dump_galaxy( std::string path)
{
	filestream.open(path, std::ios::out);
	for (size_t i = 0; i < count; i++) {
		filestream << i << ';';
		filestream << coord[i][0] << ';' << coord[i][1] << ';' << coord[i][2] << ';';
		filestream << velosites[i][0]<< ';' << velosites[i][1] << ';' << velosites[i][2] << ';';
		filestream << mass[i] << '\n';
	}
	filestream.close();
	
}

// «аписывает в файл значени€ отклонений от законов сохранени€ на текущем шаге от предыдущего без вычислени	
void NBodyData::dump_errors(std::string path)
{
	static bool first_run = true;
	std::vector<vector3> correction_vector(3, vector3());
	std::vector<value_type> correction_scalar(2, 0.);
	if (first_run) {
		first_run = false;
		std::fstream clear_file(path, std::ios::out);
		clear_file.close();
	}
	filestream.open(path, std::ios_base::app);
	filestream << prev_tottal_potential_energy - current_total_potential_energy << SEPARATOR
		<< prev_total_kinetic_energy - current_total_kinetic_energy << SEPARATOR
		<< (prev_total_kinetic_energy - current_total_kinetic_energy) +
		(prev_tottal_potential_energy - current_total_potential_energy) << std::endl;
	
	filestream.close();



	
}

vector3 NBodyData::calculate_total_impulse()
{
	prev_total_impulse = current_total_impulse;
	current_total_impulse = vector3();
	vector3 cor;
	#pragma omp parallel for
	for (int body = 0; body < count; body++) {
		current_total_impulse = Kahan_sum(current_total_impulse, velosites[body] * mass[body], &cor);
	}
	return current_total_impulse;
}

// –ассчитывает потенциальную энергию и сохран€ет предыдущее значение, возвращает новое значение
value_type NBodyData::calculate_total_potential_energy()
{
	prev_tottal_potential_energy = current_total_potential_energy;
	current_total_potential_energy = 0.;
	for (size_t body1 = 0; body1 < count; body1++) {
		value_type pot_energy = 0.;
		value_type cor = 0.;
		for (size_t body2 = 0; body2 < count; body2++) {
			if (body1 == body2)
				continue;
			current_total_potential_energy = Kahan_sum(current_total_potential_energy, potential_energy(body1, body2), &cor);
		}
		// current_total_potential_energy += Kahan_sum(current_total_potential_energy, pot_energy, &cor);
	}
	return current_total_potential_energy;
}

value_type NBodyData::calculate_total_kinetic_energy()
{
	prev_total_kinetic_energy = current_total_kinetic_energy;
	current_total_kinetic_energy = 0.;
	value_type cor = 0.;
	for (size_t body = 0; body < count; body++) {
		current_total_kinetic_energy = Kahan_sum(current_total_kinetic_energy, mass[body] * velosites[body].norm() / 2, &cor);
	}
	return current_total_kinetic_energy;
}

void NBodyData::calculate_total_energy()
{
	calculate_total_kinetic_energy();
	calculate_total_potential_energy();
}

vector3 NBodyData::calculate_total_force(vector3 d_coord, size_t id)
{
	vector3 total_force;
	vector3 correction;
	for (size_t body2 = 0; body2 < count; body2++) {
		if(id != body2)
			total_force = Kahan_sum(total_force, force(coord[id] + d_coord, coord[body2], mass[id], mass[body2]), &correction);
	}
	return total_force;
}





size_t NBodyData::get_count() const
{
	return count;
}

void NBodyData::generate_galaxy(vector3 center, value_type radius, value_type total_mass, size_t count)
{
	add_body(center, vector3(), total_mass * 0.999);
	value_type body_mass = (total_mass - mass[0]) / (count - 1);
	for (size_t body = 1; body < count; body++) {
		value_type direction = rand() % 2 == 0 ? 1. : -1.;
		vector3 norm(0, 0, direction);
		//vector3 norm(-radius, radius);
		vector3 new_body(-radius, radius);
		norm = new_body ^ norm;
		norm = norm * vector3(0.1, 0.15);
		value_type velosity_multiplicator = sqrt( G * mass[0] / (new_body - center).length() );
		vector3 velosity = norm / norm.length() * velosity_multiplicator;
		add_body(new_body, velosity, body_mass);

	}

}

vector3 NBodyData::total_impulse() const
{
	return current_total_impulse;
}

value_type NBodyData::total_energy() const
{
	return current_total_kinetic_energy + current_total_potential_energy;
}

vector3 NBodyData::last_total_impulse() const
{
	return prev_total_impulse;
}

value_type NBodyData::last_total_energy() const
{
	return prev_total_kinetic_energy + prev_tottal_potential_energy;
}

value_type NBodyData::impulce_err() const
{
	return 100 * ((last_total_impulse() - total_impulse()).length() / total_impulse().length());
}

value_type NBodyData::energy_err() const
{
	return fabs(100 * (last_total_energy() - total_energy()) / total_energy());
}

//vector3 NBodyData::total_impulse() const
//{
//	for (size_t body = 0; body < count; body++) {
//		
//	}
//	return vector3();
//}

//template<class WorkType>
//inline WorkType NBodyData::Kahan_sum(WorkType a, WorkType b, WorkType* correction)
//{
//	WorkType	corrected = b - *correction;
//	WorkType	new_sum = a + corrected;
//	*correction = (new_sum - a) - corrected;
//	return new_sum;
//}

//template<class V>
//V NBodyData::Kahan_sum_array(V a, V b)
//{
//	V correction(0.);
//	V sum = a[0];
//	for (size_t body = 0; body < count; body++) {
//		sum
//	}
//}
