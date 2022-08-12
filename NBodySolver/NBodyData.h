#pragma once

#include <vector>
#include<fstream>
#include<string>
#include "Point3.h"

// Для задачи трех тел в равнобедренном треушольнике показал хороший результат по траектории (1e-7)
#define NBODY_MIN_RADIUS 1e-5
#define G 1
#define SEPARATOR ';'

typedef double value_type;
typedef Point3<value_type> vector3;

class NBodyData
{
	size_t count;
	value_type time;
	size_t step;




	std::vector<vector3> coord;
	std::vector<vector3> velosites;
	std::vector<value_type> mass;
	// std::vector<value_type> a;

	vector3 total_impulse;
	vector3 total_impulse_moment;
	vector3 mass_center;
	value_type current_total_kinetic_energy;
	value_type prev_total_kinetic_energy;
	value_type current_total_potential_energy;
	value_type prev_tottal_potential_energy;

	std::fstream filestream;

public:
	NBodyData();

	vector3 force(const vector3& v1, const vector3& v2, const value_type mass1, const value_type mass2);
	void increase_time(value_type dt);

	value_type potential_energy(size_t vody1, size_t body2) const;
	//void add_body( const vector3& r, const vector3& v, const value_type& m, const value_type& a );
	value_type get_time() const;
	size_t get_step() const;
	vector3* get_coord();
	const vector3* get_coord() const;
	vector3* get_velosites();
	const vector3* get_velosites() const;
	const value_type* get_mass() const;

	void load_galaxy(const char* path);
	void add_body(std::vector<value_type> param);
	void add_body(vector3 coord, vector3 velosites, value_type mass);

	void dump_galaxy( std::string path);
	void dump_errors(std::string path);

	value_type calculate_total_potential_energy();
	value_type calculate_total_kinetic_energy();
	//void dump_statistics();
	//void dump_body( size_t n );
	size_t get_count() const;

	void generate_galaxy( vector3 center, value_type radius, value_type total_mass, size_t count );

	//vector3 total_impulse() const;
	//vector3 total_impulce_moment() const;
	//vector3 mass_center() const;
	value_type total_energy() const;
	//vector3 last_total_impulce() const;
	//vector3 last_total_impulce_moment() const;
	//vector3 last_mass_center() const;
	value_type last_total_energy() const;

	//value_type impulce_err() const;
	//value_type impulce_moment_err() const;
	value_type energy_err() const;
	template<class WorkType>
	WorkType Kahan_sum(WorkType a, WorkType b, WorkType* correction) {
		WorkType	corrected = b - *correction;
		WorkType	new_sum = a + corrected;
		*correction = (new_sum - a) - corrected;
		return new_sum;
	}
	//template<class V>
	// V Kahan_sum_array(V a, V b);


};

class MyCharCType : public std::ctype<char>
{
public:
	mask const* get_table()
	{
		static std::vector<std::ctype<char>::mask> table(classic_table(),
			classic_table() + table_size);
		table[';'] = space;
		return &table[0];
	}

	MyCharCType(size_t refs = 0)
		: std::ctype<char>(get_table(), false, refs)
	{
	}
};

