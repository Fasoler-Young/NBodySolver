#include<iostream>
#include "NBodeSolverEuler.h"
#include "NBodySolverRungeKutta.h"

#include <cstring>
#include <iostream>
#include <string>

#include<ctime>

#include"FileManager.h"

int main() {
	unsigned int start_time = clock();
	NBodyData data;
	FileManager fm("test.conf", &data);
	NBodySolverRungeKutta solver(&data);
	value_type total_step_time = 0.;
	value_type total_write_time = 0.;
	value_type step_time_start = 0.;
	value_type write_time_start = 0.;
	value_type write_time_end = 0.;
	value_type step_time_end = 0.;
	vector3 time_calculating;
	std::string name = "";
	data.calculate_total_energy();
	data.calculate_total_impulse();
	value_type start_energy = data.total_energy();
	value_type start_impulse = data.total_impulse().length();
	bool need_dump = false;
	while (data.get_time() < fm.get_end_time()) {


		// Подгоняем следующий шаг под запись
		if (data.get_time() + *(fm.get_dt()) > fm.next_dump_time(data.get_time())) {
			// С подгоном шага возникают трудности при работе с постоянным временем, 
			// можно было оставить это только для адаптивных методов, но пока рано. 
			// На данном этапе создандим флаг для понимания, что на следующем нужен дамп
			//fm.set_dt(fm.next_dump_time(data.get_time()) - data.get_time());
			need_dump = true;
			time_calculating.y = clock();
			data.calculate_total_energy();
			data.calculate_total_impulse();
			time_calculating.z = clock();
			time_calculating.x += time_calculating.z - time_calculating.y;
		}
		else if (need_dump) {
			need_dump = false;
			write_time_start = clock();
			fm.dump_galaxy(&data);
			write_time_end = clock();
			total_write_time += write_time_end - write_time_start;
			time_calculating.y = clock();
			data.calculate_total_energy();
			data.calculate_total_impulse();
			time_calculating.z = clock();
			time_calculating.x += time_calculating.z - time_calculating.y;

			fm.dump_errors(&data);
			std::cout << int(data.get_time() / fm.get_dump_time()) << ' ' << size_t(data.get_time()*100)/100. << '\t'
				<< "step time: " << total_step_time 
				<< " write_time: " << total_write_time 
				<< " calculate error time: " << time_calculating.x
				<< std::endl;
		} // Считаем энергию на предыдущем записи шаге, чтобы увидеть динамику

		write_time_start = clock();


		// fm.dump_errors(&data);
		
		write_time_end = clock();
		total_write_time += write_time_end - write_time_start;
		step_time_start = clock();
		solver.step(fm.get_dt());
		step_time_end = clock();
		total_step_time += step_time_end - step_time_start;

	}
	data.calculate_total_energy();
	data.calculate_total_impulse();
	std::cout << fabs(start_energy - data.total_energy()) << '\n'
		<< fabs(start_impulse - data.total_impulse().length()) << '\n'
		<< "\ncalculate err: " << time_calculating.x 
		<< "\nstep time: " << total_step_time
		<< "\nwrite time: " << total_write_time << std::endl;

	getchar();
	return 0;
}