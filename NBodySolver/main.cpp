#include<iostream>
#include "NBodeSolverEuler.h"
#include "NBodySolverRungeKutta.h"
#include "NBodySolverDormanPrince.h"

#include <cstring>
#include <iostream>
#include <string>

#include<ctime>

#include"FileManager.h"

int timer(void (*funk)(NBodyData* data), NBodyData* data);
int timer(void (*funk)());

int main() {
	unsigned int start_time = clock();
	NBodyData data;
	FileManager fm("test.conf", &data);
	NBodySolverDormanPrince solver(&data);
	std::string name = "";
	//data.calculate_total_energy();
	//data.calculate_total_impulse();
	//fm.calculate_condition(&data);
	value_type start_energy = data.total_energy();
	value_type start_impulse = data.total_impulse().length();
	bool need_dump = false;
	fm.dump_galaxy(&data);
	while (data.get_time() < fm.get_end_time()) {

		// Подгоняем следующий шаг под запись
		if (data.get_time() + *(fm.get_dt()) > fm.next_dump_time(data.get_time())) {
			// С подгоном шага возникают трудности при работе с постоянным временем, 
			// можно было оставить это только для адаптивных методов, но пока рано. 
			// На данном этапе создандим флаг для понимания, что на следующем нужен дамп
			//fm.set_dt(fm.next_dump_time(data.get_time()) - data.get_time());
			need_dump = true;
			fm.calculate_condition(&data);
		}
		else if (need_dump) {
			need_dump = false;
			fm.calculate_condition(&data);
			fm.dump_galaxy(&data);
			fm.dump_errors(&data);
			fm.log(&data, start_energy, start_impulse);
		} // Считаем энергию на предыдущем записи шаге, чтобы увидеть динамику
		
		fm.step();
		//fm.calculate_condition(&data);

	}
	fm.dump_galaxy(&data);
	fm.dump_errors(&data);
	fm.log(&data, start_energy, start_impulse);
	//data.calculate_total_energy();
	//data.calculate_total_impulse();
	//std::cout << fabs(start_energy - data.total_energy()) << '\n'
	//	<< fabs(start_impulse - data.total_impulse().length()) 
	//	<< "\nsteps: " << data.get_step()
	//	<< "\nmethod: " << fm.get_method_name()
	//	<< std::endl;

	getchar();
	return 0;
}

int timer(void(*funk)(NBodyData* data), NBodyData* data)
{
	int start = clock();
	funk(data);
	return clock() - start;
}

int timer(void(*funk)())
{
	int start = clock();
	funk();
	return clock() - start;
}
