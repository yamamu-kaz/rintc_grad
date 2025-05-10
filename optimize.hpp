/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef main_h
#define main_h

#include <string>

namespace rintdwr 
{
std::string get_current_datetime();

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
	const std::string data_file);

template<typename RealScalar>
void calc_bppm(
	const std::string param_file,
	const std::string data_file,
	const std::string outputfile="");

}
#endif /* main_h */
