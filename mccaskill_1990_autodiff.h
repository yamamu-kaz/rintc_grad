/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_MCCASKILL_1990_AUTODIFF_H_
#define RINTDWR_MCCASKILL_1990_AUTODIFF_H_

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
#include <cctype>
#include <utility>

#include"interval_type.h"
#include"real_logsumexp.h"

#include "optimize.hpp"

namespace rintdwr {

typedef double DiffFloat;
typedef std::pair<bool, DiffFloat>MemDiff;


template<typename T>
void print2DVector(const std::vector<std::vector<T>>& vec);

template<typename T>
void print1DVector(const std::vector<T>& vec);

template<typename RealScalar>
class SimpleMcCaskillAutoDiff {

public:
	typedef std::pair<bool, RealScalar> MemReal;
	typedef std::vector<std::vector<RealScalar>> BppmType;

	SimpleMcCaskillAutoDiff(
		// const std::string sequence,
		const std::string param_file_name,
		const double temperature,
		// const int max_span,
		const int max_loop
	);

	std::pair<std::vector<std::vector<RealScalar>>, std::vector<std::vector<bool>>> run();
	std::vector<std::vector<DiffFloat>> diff(const std::vector<std::vector<bool>>& isbppm_clipped);

	void train(const int nepoch, 
			   const double learning_rate, 
			   int nsample, 
			   const int sequence_length,
			   int max_span,
			   const std::string logfile,
			   const double logstack_init,
			   const double weight_central,
			   const std::string datafile);
	
	void calc(
				const std::string datafile,
				const std::string outputfile=""
	);

	void setSequence(const std::string sequence,const int max_span);
	void init6Vectors();
	void init6DiffVectors();
	void reset_dlogstack();
	void add_dlogstack(const DiffFloat value, const int i, const int j);
	void print_dlogstack();
	std::vector<RealScalar> calc_bppmsum(
		const std::vector<std::vector<RealScalar>>& bppm,
		const std::string seq,
		const int max_span);
	std::vector<DiffFloat>calc_dbppmsum(
		const std::vector<std::vector<DiffFloat>>& dbppm_dp,
		const std::string seq,
		const int max_span
	);
	void setPidx(const int pidx1, const int pidx2);
private:
	std::string sequence;
	std::string param_file_name;
	double temperature;
	int max_span;
	int max_loop;
	int n;

	std::vector<MemReal>Z; 
	std::vector<MemReal>Z1;
	std::vector<MemReal>Zb;
	std::vector<MemReal>Zm;
	std::vector<MemReal>Zm1;
	std::vector<MemReal>Wb;

	std::vector<MemDiff>dZ_dp; 
	std::vector<MemDiff>dZ1_dp;
	std::vector<MemDiff>dZb_dp;
	std::vector<MemDiff>dZm_dp;
	std::vector<MemDiff>dZm1_dp;
	std::vector<MemDiff>dWb_dp;

	int at(const int i, const int j);

	RealScalar GetZ(const int i, const int j);
	RealScalar GetZ1(const int i, const int j);
	RealScalar GetZb(const int i, const int j);
	RealScalar GetZm(const int i, const int j);
	RealScalar GetZm1(const int i, const int j);
	RealScalar GetWb(const int i, const int j);

	const int pidx1_lower = 1;
	const int pidx2_lower = 1;
	int pidx1 = 1;
	int pidx2 = 1;
	static constexpr int pidx1Size = 12;
	static constexpr int pidx2Size = 12;
	DiffFloat dlogstack[pidx1Size][pidx2Size] = {0.0};
	bool update_mask[pidx1Size][pidx2Size] = {0};
	

	DiffFloat Get_dZ_dp(const int i, const int j);
	DiffFloat Get_dZ1_dp(const int i, const int j);
	DiffFloat Get_dZb_dp(const int i, const int j);
	DiffFloat Get_dZm_dp(const int i, const int j);
	DiffFloat Get_dZm1_dp(const int i, const int j);
	DiffFloat Get_dWb_dp(const int i, const int j);

	void initVector(std::vector<MemReal>&vec, const int length, const MemReal value);
	void write_bppm(
		const std::string file,
		const std::vector<std::vector<RealScalar>> &bppm,
		const std::string seq);

	int load_update_mask(const std::string file);
	std::vector<std::pair<std::string, BppmType>> load_sequence_and_bppm(const std::string datafile);
};

}

#endif//RINTDWR_MCCASKILL_1990_AUTODIFF_H_