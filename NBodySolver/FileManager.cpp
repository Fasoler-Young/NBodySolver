#include "FileManager.h"

FileManager::FileManager(std::string path_conf, NBodyData* data)
{
	
	// Читаем конфигурацию для определения путей
	paths["path_config"] = path_conf;
	filestream.open(paths["path_config"], std::ios::in);
	std::string key;
	std::string value;
	while (filestream >> key) {
		filestream >> paths[key];
	}
	filestream.close();

	
	// Флаг о том, что задача генерируется случайно
	if(paths["path_data_input"] == "generate"){
		dt = 1e-5;
		end_time = 2.5;
		output_files_count = 151;
		data->generate_galaxy(vector3(), 20., 10., 50);
	}
	else {
		// Читаем данные задачи
		auto loc = std::locale(std::locale(), new MyCharCType());
		filestream.imbue(loc);

		filestream.ignore(std::numeric_limits<std::streamsize>::max(),
			filestream.widen('\n'));
		filestream.open(paths["path_data_input"], std::ios::in);
		value_type data1;

		// Начальный шаг и время расчетов индивидуальны для задачи
		filestream >> dt >> end_time >> output_files_count;

		std::vector<value_type> param;
		while (filestream >> data1)
		{
			param.push_back(data1);
			if (param.size() == 7) {
				data->add_body(param);
				param.clear();
			}
		}
		filestream.close();
	}
	dump_step = ceil(end_time / dt / output_files_count);
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

value_type FileManager::get_output_files_count() const
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

void FileManager::dump_galaxy(NBodyData* data)
{
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
}

void FileManager::dump_errors(NBodyData* data)
{
	static bool first_run = true;
	std::vector<vector3> correction_vector(3, vector3());
	std::vector<value_type> correction_scalar(2, 0.);
	if (first_run) {
		first_run = false;
		std::fstream clear_file(get_path_data_error(), std::ios::out);
		clear_file.close();
	}
	filestream.open(get_path_data_error(), std::ios_base::app);
	filestream /*<< prev_tottal_potential_energy - current_total_potential_energy << SEPARATOR
		<< prev_total_kinetic_energy - current_total_kinetic_energy << SEPARATOR
		<< (prev_total_kinetic_energy - current_total_kinetic_energy) +
		(prev_tottal_potential_energy - current_total_potential_energy) << std::endl;*/
		<< data->energy_err() << SEPARATOR
		<< data->impulce_err()
		<< std::endl;

	filestream.close();
}

void FileManager::set_dt(value_type new_t)
{
	dt = new_t;
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
