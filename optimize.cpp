/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

#include"misc.h"
#include"rintx.h"
#include"mccaskill_1990_autodiff.h"
#include"centroid_fold.h"


#ifdef _WIN64 
#include<windows.h>
#endif

namespace rintdwr {

std::string get_current_datetime() {
  // Get current time in seconds since epoch
  time_t rawtime;
  time(&rawtime);

  // Convert time to tm struct for manipulation
  struct tm * timeinfo = localtime(&rawtime);

  // Create string stream for formatted output
  std::stringstream ss;

  // Use setw to ensure two-digit formatting for all fields
  ss << std::setfill('0') 
     << std::setw(4) << timeinfo->tm_year + 1900
     << std::setw(2) << timeinfo->tm_mon + 1
     << std::setw(2) << timeinfo->tm_mday << "_"
     << std::setw(2) << timeinfo->tm_hour
     << std::setw(2) << timeinfo->tm_min
     << std::setw(2) << timeinfo->tm_sec;

  // Extract the formatted string from the stream
  return ss.str();
}

// Explicit instantiation for double
template void optimize<Floating>(
	const int max_loop, 
	const int max_span, 
	const int nsample, 
	const int sequence_length, 
	const int nepoch, 
	const double learning_rate,
	const double logstack_init,
	const double weight_central,
	const std::string param_file,
	const std::string data_file) ;

template<typename RealScalar>
void optimize(
	const int max_loop, 
	const int max_span, 
	const int nsample, 
	const int sequence_length, 
	const int nepoch, 
	const double learning_rate,
	const double logstack_init,
	const double weight_central,
	const std::string param_file,
	const std::string data_file) 
{

	//random seq
	std::cout << "\n<-- start energy parameter optimization experiment --> \n";

	const double temperature = 37.0;
	// const std::string param_file_name = "modified_nuculeobases";
	// const int max_loop = 30;
	// const int max_span = 20;
	// const int nsample = 10;
	// const int sequence_length = 200;
	// const int nepoch = 20;
	// const double learning_rate = 1.0;
	const std::string logfile = "./logs/rintc-opt." + get_current_datetime() + ".txt";

	std::ofstream outfile(logfile, std::ios::app);
	if (outfile.is_open()) {
		outfile << "rintc-opt stacking parameter optimization\n";
		outfile << "temperature            = " << temperature << "\n"; 
		// outfile << "param_file_name        = " << param_file_name << "\n";
		outfile << "max_loop               = " << max_loop << "\n";
		outfile << "max_span               = " << max_span << "\n";
		outfile << "nsample                = " << nsample << "\n";
		outfile << "sequence_length        = " << sequence_length << "\n";
		outfile << "nepoch                 = " << nepoch << "\n";
		outfile << "learning_rate          = " << learning_rate << "\n";
		outfile << "weight_centeral        = " << weight_central << "\n";
		outfile << "logstack_init          = " << logstack_init << "\n";
		outfile << "param_file             = " << param_file << "\n";
		outfile << "data_file              = " << data_file << "\n\n";
	}

	std::string datafile = data_file=="random" ? "" : data_file;

	outfile.close();

	auto smc = SimpleMcCaskillAutoDiff<RealScalar>(param_file, temperature, max_loop);
	smc.train(nepoch, learning_rate, nsample, sequence_length, max_span, logfile, logstack_init, weight_central, datafile);

}


// Explicit instantiation for double
template void calc_bppm<Floating>(
	const std::string param_file,
	const std::string data_file,
	const std::string outputfile
	)  ;


template<typename RealScalar>
void calc_bppm(
	const std::string param_file,
	const std::string data_file,
	const std::string outputfile
	) 
{
	const double temperature = 37.0;
	const int max_loop = 30;
	auto smc = SimpleMcCaskillAutoDiff<RealScalar>(param_file,temperature,max_loop );
	if(outputfile.empty()){
		smc.calc(data_file);
	}else{
		smc.calc(data_file,outputfile);
	}
}


}
