#pragma once
#include<string>
#include<map>
#include<iostream>
#include<ctime>

#include"NBodyData.h"

#include"NBodySolver.h"
#include"NBodeSolverEuler.h"
#include"NBodySolverRungeKutta.h"
#include"NBodySolverDormanPrince.h"
#include"NBodySolverAdamsBashfort.h"


class FileManager
{
private:
	std::map<std::string, std::string> paths;
	value_type dt;
	value_type end_time;
	size_t dump_step;
	value_type dump_time;
	size_t output_files_count;
	NBodySolver* solver;
	


	size_t total_step_time, total_write_time, total_calculate_err_time;


	std::fstream filestream;

	void create_solvers(NBodyData* data);
	void create_galaxy(NBodyData* data);

public:
	FileManager(std::string path_conf, NBodyData* data);

	value_type get_end_time() const;
	size_t get_dump_step() const;
	value_type get_dump_time() const;
	size_t get_output_files_count() const;
	value_type* get_dt();
	value_type next_dump_time(value_type cur_time);
	
	void step();
	void calculate_condition(NBodyData* data);

	void dump_galaxy(NBodyData* data);
	void dump_errors(NBodyData* data);
	void log(NBodyData* data, value_type E_0, value_type P_0);

	void set_dt(value_type new_t);
	std::string get_method_name();
	std::string get_path_log();
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

