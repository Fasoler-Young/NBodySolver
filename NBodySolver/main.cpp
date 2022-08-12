#include<iostream>
#include "NBodeSolverEuler.h"

#include <cstring>
#include <iostream>
#include <string>

#include<ctime>

#include"FileManager.h"

int main() {
	unsigned int start_time = clock();
	NBodyData data;
	FileManager fm("test.conf", &data);
	//const char* path = "test.csv";
	//const char* path = "triangle_x_y.dat";
	//const char* path = "two_equals_body_circle.dat";
	//data.load_galaxy(path);
	//unsigned int end_time = clock();
	//value_type dt = 1e-6;
	NBodeSolverEuler solver(&data);
	//std::string path2 = "G:\\NBody\\out\\out_";
	//std::string path_err = "G:\\NBody\\err\\err_";
	//std::string end = ".dat";
	//value_type stop = 100.;
	//size_t dump_step = ceil( stop/dt/151);
	value_type total_step_time = 0.;
	value_type total_write_time = 0.;
	value_type step_time_start = 0.;
	value_type write_time_start = 0.;
	value_type write_time_end = 0.;
	value_type step_time_end = 0.;
	std::string name = "";

	while (data.get_time() < fm.get_end_time()) {
		if (data.get_step() % fm.get_dump_step() == 0) {
			write_time_start = clock();
			
			//name = (path2 + (std::to_string(data.get_step() / dump_step)) + end);
			fm.dump_galaxy(&data);
			//data.dump_galaxy(fm.get_path_data_output());
			write_time_end = clock();
			total_write_time += write_time_end - write_time_start;
			//value_type pot_en = -G / (data.get_coord()[1] - data.get_coord()[0]).length();
			//value_type kin_en = data.get_velosites()[1].norm() / 2;
			data.calculate_total_kinetic_energy();
			data.calculate_total_potential_energy();
			fm.dump_errors(&data);
			//data.dump_errors(path_err + (std::to_string(0)) + end);
			std::cout << data.get_step() / fm.get_dump_step() << ' ' << size_t(data.get_time()) << '\t'
				//<< "x: " << data.get_coord()[1][0] << "\t"
				//<< "y: " << data.get_coord()[1][1] << "\t"
				//<< "z: " << data.get_coord()[1][2] << "\t"
				/*<< "v: " << sqrt(data.get_velosites()[1].norm()) << "\t"
				<< "pot_en: " << pot_en << "\t"
				<< "kin_in: " << kin_en << "\t"
				<< "en: " << pot_en + kin_en*/
				<< data.energy_err()
				<< std::endl;
		}

		step_time_start = clock();
		solver.step(fm.get_dt());
		step_time_end = clock();
		total_step_time += step_time_end - step_time_start;

	}
	std::cout << "step time: " << total_step_time
		<< "\nwrite time: " << total_write_time << std::endl;

	getchar();
	return 0;
}