/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_PARAMETER_H_
#define RINTDWR_PARAMETER_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <functional>
#include <complex>
#include <random>
#include <chrono>
#include <sstream>

namespace rintdwr {

namespace parasor_param {

int GetBaseType(const char c);
int GetPairType(const char c1, const char c2);
int GetPairTypeReverse(const char c1, const char c2);

double ParMultiloopClosing();
double ParMultiloopInternal();
double ParDangling(const int type, const int five, const int three, const bool ext, const std::string& sequence);
double ParHairpinEnergy(const int i, const int j, const std::string& sequence);
double ParLoopEnergy(const int i, const int j, const int p, const int q, const std::string& sequence);
void InitializeParameter(const std::string& parameter_file_name, const double temperature);

double dparstack_dp(const int i, const int j);
double dParLoopEnergy_dp(
	const int i, 
	const int j, 
	const int p, 
	const int q, 
	const std::string& sequence, 
	const int pidx1, 
	const int pidx2);
void update_logstack(const double dif, const int i, const int j);
void update_parstack();

template <size_t pidx1Size, size_t pidx2Size>
void reset_logstack(const double value, const double weight_central, bool (&mask)[pidx1Size][pidx2Size]);

template <size_t pidx1Size, size_t pidx2Size>
void print_logstack(const std::string logfile,  bool (&mask)[pidx1Size][pidx2Size]);


extern bool counting;

}

}


#endif//RINTDWR_PARAMETER_H_