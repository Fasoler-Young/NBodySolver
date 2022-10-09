#include "FileManager.h"

void FileManager::create_solvers(NBodyData* data)
{
	std::string solver_name;
	filestream >> solver_name;
	// Нет смылсла оптимизировать, компилятор все равно удаляет ненужные
	NBodySolverEuler* euler = new NBodySolverEuler(data);
	RK4* runge_kutta = new RK4(data);
	NBodySolverDormanPrince* dormane_prince = new NBodySolverDormanPrince(data);
	NBodySolverAdamsBashfort* adams = new NBodySolverAdamsBashfort(data);
	if (!solver_name.compare(euler->method_name()))
		solver = euler;
	else if (!solver_name.compare(runge_kutta->method_name())) {
		value_type local_err;
		filestream >> local_err;
		runge_kutta->set_max_local_err(local_err);
		solver = runge_kutta;
	}
	else if (!solver_name.compare(dormane_prince->method_name())) {
		value_type local_err;
		filestream >> local_err;
		dormane_prince->set_max_local_err(local_err);
		solver = dormane_prince;
	}
	else if (!solver_name.compare(adams->method_name())) {
		solver = adams;
	}
	std::cout << "method from file: " << solver_name << "\t current method: " << solver->method_name() << std::endl;

		
}

// !!! Метод пока неустойчив к ошибкам в файле и надо это иметь в виду
void FileManager::create_galaxy(NBodyData* data)
{
	std::string type;
	filestream >> type;
	if (!type.compare("generate")) {
		vector3 center;
		value_type radius, total_mass;
		size_t count;
		filestream >> center.x >> center.y >> center.z
			>> radius >> total_mass >> count;
		value_type v = 0;// -radius * 1.75; //2.23606797749979;
		data->generate_galaxy(center/* - vector3(1.6*radius, 0, 0)*/, radius, total_mass, count, vector3(0, -v, 0));
		//data->generate_galaxy(center + vector3(1.6 * radius, 0, 0), radius, total_mass, count, vector3(0, v, 0));
	}
	else {
		std::vector<value_type> param;
		value_type tmp;
		while (filestream >> tmp)
		{
			param.push_back(tmp);
			if (param.size() == 7) {
				data->add_body(param);
				param.clear();
			}
		}
	}

}

FileManager::FileManager(std::string path_conf, NBodyData* data)
{
	total_calculate_err_time = total_step_time = total_write_time = 0;


	// Читаем конфигурацию для определения путей
	paths["path_config"] = path_conf;
	filestream.open(paths["path_config"], std::ios::in);
	std::string key;
	std::string value;
	while (filestream >> key) {
		filestream >> paths[key];
	}

	filestream.close();

	
	//// Флаг о том, что задача генерируется случайно
	//if(paths["path_data_input"] == "generate"){
	//	dt = 1e-8;
	//	end_time = 2.5;
	//	output_files_count = 151;
	//	data->generate_galaxy(vector3(), 20., 10., 10);
	//}
	//else {
	// 
		// Читаем данные задачи
	auto loc = std::locale(std::locale(), new MyCharCType());
	filestream.imbue(loc);

	filestream.ignore(std::numeric_limits<std::streamsize>::max(),
		filestream.widen('\n'));
	filestream.open(paths["path_data_input"], std::ios::in);

	// Начальный шаг и время расчетов индивидуальны для задачи
	filestream >> dt >> end_time >> output_files_count;
	create_solvers(data);
	create_galaxy(data);


	filestream.close();
	calculate_condition(data);
	E_0 = data->total_energy();
	P_0 = data->total_impulse().length();
	//}
	//char tmp[255];
	//sprintf_s(tmp, "%.e", dt);
	//paths["path_data_output"] = get_path_data_output()
	//	.replace(get_path_data_output().find("%dt"), 3, tmp)
	//	.replace(get_path_data_output().find("%method"), 7, solver->method_name());
	//paths["path_data_error"] = get_path_data_error()
	//	.replace(get_path_data_error().find("%dt"), 3, tmp)
	//	.replace(get_path_data_error().find("%method"), 7, solver->method_name());
	//paths["path_log"] = get_path_log()
	//	.replace(get_path_log().find("%dt"), 3, tmp)
	//	.replace(get_path_log().find("%method"), 7, solver->method_name());

	dump_step = size_t(ceil(end_time / dt / output_files_count));
	dump_time = end_time / output_files_count;

}

value_type FileManager::get_end_time() const
{
	return end_time;
}

size_t FileManager::get_dump_step() const
{
	return dump_step;
}

value_type FileManager::get_dump_time() const
{
	return dump_time;
}

size_t FileManager::get_output_files_count() const
{
	return output_files_count;
}

value_type* FileManager::get_dt() 
{
	return &dt;
}

value_type FileManager::next_dump_time(value_type cur_time)
{
	return (int(cur_time / dump_time) + 1) * dump_time;
}

void FileManager::step()
{
	size_t start = clock();
	solver->step(&dt);
	total_step_time += clock() - start;
}

void FileManager::calculate_condition(NBodyData* data)
{
	int start = clock();
	data->calculate_total_energy();
	data->calculate_total_impulse();
	total_calculate_err_time += clock() - start;
}

void FileManager::dump_galaxy(NBodyData* data)
{
	size_t start = clock();
	static size_t file_number = 0;
	
	filestream.open(get_path_data_output().replace(get_path_data_output().find("%i"), 2, std::to_string(file_number)), std::ios::out);
	for (size_t i = 0; i < data->get_count(); i++) {
		filestream << i << ';';
		filestream << data->get_coord()[i][0] << ';' << data->get_coord()[i][1] << ';' << data->get_coord()[i][2] << ';';
		filestream << data->get_velosites()[i][0] << ';' << data->get_velosites()[i][1] << ';' << data->get_velosites()[i][2] << ';';
		filestream << data->get_mass()[i] << '\n';
	}
	file_number++;
	filestream.close();
	total_write_time += clock() - start;
}

void FileManager::dump_errors(NBodyData* data)
{
	size_t start = clock();
	static bool first_run = true;
	if (first_run) {
		first_run = false;
		std::fstream clear_file(get_path_data_error(), std::ios::out);
		clear_file.close();
	}
	filestream.open(get_path_data_error(), std::ios_base::app);
	filestream 
		<< data->energy_err() << SEPARATOR							// 1
		<< data->impulce_err() << SEPARATOR							// 2
		<< data->total_energy() << SEPARATOR						// 3
		<< data->total_impulse().length() << SEPARATOR				// 4
		<< E_0 - data->total_energy() << SEPARATOR					// 5
		<< P_0 - data->total_impulse().length() << SEPARATOR		// 6
		<< data->get_potential_energy() << SEPARATOR				// 7
		<< data->get_kinetik_energy() << SEPARATOR					// 8
		<< fabs((E_0 - data->total_energy()) / E_0) << SEPARATOR	// 9

		<< std::endl;

	filestream.close();
	total_write_time += clock() - start;
}

void FileManager::log(NBodyData* data, value_type E_0, value_type P_0)
{
	size_t start = clock();
	static bool first_run = true;
	bool last_log = data->get_time() > end_time;
	if (first_run) {
		first_run = false;
		std::fstream clear_file(get_path_log(), std::ios::out);
		clear_file << "dump\t" << "time\t" << "E_0-E_cur\t" << "P_0-P_cur\t" << "tot_stp_time\t" << "tot_wr_time\t" << "tot_clc_err_time\n";
		std::cout  << "dump\t" << "time\t" << "E_0-E_cur\t" << "P_0-P_cur\t" << "tot_stp_time\t" << "tot_wr_time\t" << "tot_clc_err_time\n";
		clear_file.close();
	}
	filestream.open(get_path_log(), std::ios_base::app);
	filestream << size_t(data->get_time() / dump_time) << '\t'
		<< size_t(data->get_time() * 100) / 100. << '\t'
		<< fabs(E_0 - data->total_energy()) << '\t'
		<< fabs(P_0 - data->total_impulse().length()) << "\t \t"
		<< total_step_time << "\t\t"
		<< total_write_time << "\t\t"
		<< total_calculate_err_time << '\t'
		<< dt
		<< std::endl;
	if (last_log) {
		filestream << fabs(E_0 - data->total_energy()) << '\n'
			<< fabs(P_0 - data->total_impulse().length()) << '\n'
			<< data->get_number_of_function_calls() << '\n'
			<< data->get_step() << '\n'
			<< solver->method_name();

	}
	filestream.close();
	std::cout << size_t(data->get_time() / dump_time) << '\t'
		<< size_t(data->get_time() * 100) / 100. << '\t'
		<< fabs(E_0 - data->total_energy()) << '\t'
		<< fabs(P_0 - data->total_impulse().length()) << "\t \t"
		<< total_step_time << "\t\t"
		<< total_write_time << "\t\t"
		<< total_calculate_err_time << '\t'
		<< dt << std::endl;
	if (last_log) {
		std::cout << fabs(E_0 - data->total_energy()) << '\n'
			<< fabs(P_0 - data->total_impulse().length()) << '\n'
			<< data->get_number_of_function_calls() << '\n'
			<< data->get_step() << '\n'
			<< solver->method_name();

	}
	total_write_time += clock() - start;

}

void FileManager::set_dt(value_type new_t)
{
	dt = new_t;
}

std::string FileManager::get_method_name()
{
	std::string tmp = solver->method_name();
	return tmp;
}

std::string FileManager::get_path_log()
{
	return paths.at("path_log");
}

std::string FileManager::get_path_config() const
{
	return paths.at("path_config");
}

std::string FileManager::get_path_data_output()
{
	return paths.at("path_data_output");
}

std::string FileManager::get_path_data_error()
{
	return paths.at("path_data_error");
}
