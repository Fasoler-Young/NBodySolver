#pragma once
#include<string>
#include<map>
#include<fstream>
#include"NBodyData.h"

class FileManager
{
	std::map<std::string, std::string> paths;
	value_type dt;
	value_type end_time;
	size_t dump_step;
	size_t output_files_count;

	std::fstream filestream;

public:
	FileManager(std::string path_conf, NBodyData* data);

	value_type get_end_time() const;
	size_t get_dump_step() const;
	value_type get_output_files_count() const;
	value_type get_dt() const;

	void dump_galaxy(NBodyData* data);
	void dump_errors(NBodyData* data);

	std::string get_path_config() const;
	std::string get_path_data_input() const;
	std::string get_path_data_output();
	std::string get_path_data_error();

};

//class MyCharCType : public std::ctype<char>
//{
//public:
//	mask const* get_table()
//	{
//		static std::vector<std::ctype<char>::mask> table(classic_table(),
//			classic_table() + table_size);
//		table[';'] = space;
//		return &table[0];
//	}
//
//	MyCharCType(size_t refs = 0)
//		: std::ctype<char>(get_table(), false, refs)
//	{
//	}
//};

