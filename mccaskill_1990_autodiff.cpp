/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/
#include <chrono>
#include <fstream>
#include <sstream>

#include "mccaskill_1990_autodiff.h"

#include "interval_type.h"
#include "real_logsumexp.h"
#include "parameter.h"
#include "misc.h"

namespace rintdwr
{

	bool endsWith(const std::string& mainStr, const std::string& toMatch) {
		if (mainStr.size() >= toMatch.size() && 
			mainStr.compare(mainStr.size() - toMatch.size(), toMatch.size(), toMatch) == 0)
			return true;
		else
			return false;
	}

	// Function to trim whitespace from the start of a string
	std::string ltrim(const std::string &s) {
		size_t start = s.find_first_not_of(" \t\n\r\f\v");
		return (start == std::string::npos) ? "" : s.substr(start);
	}

	// Function to trim whitespace from the end of a string
	std::string rtrim(const std::string &s) {
		size_t end = s.find_last_not_of(" \t\n\r\f\v");
		return (end == std::string::npos) ? "" : s.substr(0, end + 1);
	}

	// Function to trim whitespace from both ends of a string
	std::string trim(const std::string &s) {
		return rtrim(ltrim(s));
	}

	std::string replace_string(const std::string& str, const std::string from, const std::string to) {

		std::string out = str;
		// Find the starting position of the substring to be replaced
		size_t startPos = out.find(from);
		if (startPos != std::string::npos) { // Check if the substring is found
			out.replace(startPos, from.length(), to);
		}

		return out;
	}

	int load_m6a_dataset(const std::string datafile, std::vector<std::string>& sequences) {
		std::ifstream m6afile(datafile);
		if (!m6afile.is_open()) {
			std::cerr << "Unable to open file : " << datafile << std::endl;
			return 0;
		}

		std::string line;
		std::string seq("");
		while (std::getline(m6afile, line)) {
			if (line[0] == '>') {
				continue;
				// sequences.push_back(trim(line).substr(1));  // Start a new sequence
			}
			else {
				sequences.push_back(trim(line));
				// continue;
				// sequences.back() += trim(line);  // Append to the current sequence
			}
		}
		m6afile.close();
		return sequences.size();
	}


	// Explicit instantiation for double
	template class SimpleMcCaskillAutoDiff<Floating>;

	template <typename T>
	void print2DVector(const std::vector<std::vector<T>> &vec)
	{
		std::cout << std::scientific << std::setprecision(1);
		// Loop through each row
		for (const auto &row : vec)
		{

			// Loop through each element in the row
			for (const auto element : row)
			{
				// Print the element with a space
				std::cout << element << ", ";
			}
			// Move to a new line after each row
			std::cout << std::endl;
		}
		std::cout << std::defaultfloat;
	}

	template <typename T>
	void print1DVector(const std::vector<T> &vec)
	{
		std::cout << std::scientific << std::setprecision(1);
		// Loop through each element in the row
		for (const auto element : vec)
		{
			// Print the element with a space
			std::cout << element << ", ";
		}
		std::cout << std::endl;
		std::cout << std::defaultfloat;
	}

	template <typename RealScalar>
	SimpleMcCaskillAutoDiff<RealScalar>::SimpleMcCaskillAutoDiff(
		// const std::string sequence,
		const std::string param_file_name,
		const double temperature,
		// const int max_span,
		const int max_loop) : // sequence(sequence),
							  param_file_name(param_file_name),
							  temperature(temperature),
							  // max_span(max_span),
							  max_loop(max_loop)
	{
		parasor_param::InitializeParameter(param_file_name, temperature);

		// n = int(sequence.size());
		// init6Vectors();
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::init6Vectors()
	{
		std::cout << "initialized Z,Z1,Zb,Zm,Zm1,Wb vectors to false,0.0\n";
		this->initVector(Z, (n + 1) * 2, std::make_pair(false, RealScalar(0.0)));
		this->initVector(Z1, (n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
		this->initVector(Zb, (n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
		this->initVector(Zm, (n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
		this->initVector(Zm1, (n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
		this->initVector(Wb, (n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));

		return;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::reset_dlogstack()
	{
		std::fill(&(dlogstack[0][0]), &(dlogstack[0][0]) + this->pidx1Size*this->pidx2Size, DiffFloat(0.0));
		return;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::add_dlogstack(const DiffFloat value, const int i, const int j)
	{
		this->dlogstack[i][j] += value;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::print_dlogstack()
	{
		std::cout << "[dlogstack]\n";
		std::cout << std::scientific << std::setprecision(1);

		for (int i = 0; i < this->pidx1Size; ++i)
		{
			for (int j = 0; j < this->pidx2Size; ++j)
			{
				std::cout << this->dlogstack[i][j] << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << std::defaultfloat;
	}

	template <typename RealScalar>
	int SimpleMcCaskillAutoDiff<RealScalar>::load_update_mask(const std::string file)
	{
		std::ifstream infile(file);
		if (!infile.is_open()) {
			std::cerr << "Error opening file : " << file << std::endl;
			return 1;
		}

		std::string line;
		int row = 0;

		while (std::getline(infile, line)) {
			line = trim(line);

			if (line.find("/*") == 0) {
				// Skip the comment line
				continue;
			}

			// Parse the line to extract integers
			std::vector<std::string> result;
			std::istringstream iss(line);
			std::string token;

			// Use a string stream to extract tokens separated by whitespace
			while (iss >> token) {
				result.push_back(token);
			}

			for (int col=0; col < this->pidx2Size; ++col) {
				if (col>=result.size()) break;
				this->update_mask[row][col] = result[col][0]=='0'? 0: 1;
			}

			row++;
			if (row >= this->pidx1Size) {
				break;
			}
		}

		infile.close();

		return 0;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::write_bppm(
		const std::string file,
		const std::vector<std::vector<RealScalar>> &bppm,
		const std::string seq)
	{
		std::ofstream outfs(file, std::ios::app);
		if (!outfs.is_open()) {
			std::cerr << "cannot write to : " << file << std::endl;
			return;
		}

		outfs << seq << std::endl;
		for (size_t k=1; k<bppm.size(); ++k){
			auto& vec = bppm[k];
			for (size_t i=1; i < vec.size()-1; i++) {
				outfs << std::setw(8) << std::scientific << std::setprecision(2) << vec[i] << ", ";
			}
			outfs << vec.back() << std::endl;
    	}
		outfs << std::endl;

		outfs.close();
		return;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::init6DiffVectors()
	{
		// std::cout << "initialized dZ_dp,dZ1_dp,dZb_dp,dZm_dp,dZm1_dp,dWb_dp vectors to false,0.0\n";
		this->initVector(dZ_dp, (n + 1) * 2, std::make_pair(false, DiffFloat(0.0)));
		this->initVector(dZ1_dp, (n + 1) * (max_span + 1), std::make_pair(false, DiffFloat(0.0)));
		this->initVector(dZb_dp, (n + 1) * (max_span + 1), std::make_pair(false, DiffFloat(0.0)));
		this->initVector(dZm_dp, (n + 1) * (max_span + 1), std::make_pair(false, DiffFloat(0.0)));
		this->initVector(dZm1_dp, (n + 1) * (max_span + 1), std::make_pair(false, DiffFloat(0.0)));
		this->initVector(dWb_dp, (n + 1) * (max_span + 1), std::make_pair(false, DiffFloat(0.0)));

		return;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::initVector(std::vector<MemReal> &vec, const int length, MemReal value)
	{
		vec.resize(length);

		// Initialize each element with std::make_pair(false, RealScalar(0.0))
		for (auto &elem : vec)
		{
			elem = value; // Adjust RealScalar(0.0) as needed
		}
	}

	template <typename RealScalar>
	int SimpleMcCaskillAutoDiff<RealScalar>::at(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	}

	template <typename RealScalar>
	RealScalar SimpleMcCaskillAutoDiff<RealScalar>::GetZ(const int i, const int j)
	{
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)
			return RealScalar(1.0);
		if (i == j + 1)
			return RealScalar(1.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[index].first)
			return Z[index].second;

		// RealScalar ans = RealScalar(1.0);
		// for (int k = i; k <= j; ++k)ans += GetZ(i, k - 1) * GetZ1(k, j);
		// return (Z[i][j] = std::make_pair(true, ans)).second;

		RealScalar ans = RealScalar(0.0);
		if (i == 1)
		{

			// no base pair
			ans += RealScalar(1.0);

			// there is only one outermost base pair (i, *)
			ans += GetZ1(i, j);

			// there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
			for (int k = i; k <= j - 1; ++k)
			{
				ans += GetZ(i, k) * GetZ1(k + 1, j);
			}
		}
		else
		{
			// i is not paired
			ans += GetZ(i + 1, j);

			//(i,k) form a base pair
			for (int k = i + TURN + 1; k <= j && k - i <= max_span; ++k)
			{
				const int type = parasor_param::GetPairType(sequence[i - 1], sequence[k - 1]);
				if (type == 0)
					continue;
				ans += GetZb(i, k) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence))) * GetZ(k + 1, j);
			}
		}
		return (Z[index] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	RealScalar SimpleMcCaskillAutoDiff<RealScalar>::GetZ1(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);
		if (i == j)
			return RealScalar(0.0);

		if ((j - i) > max_span)
		{
			return GetZ1(i, i + max_span);
		}

		if (Z1[at(i, j)].first)
			return Z1[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0)
		{
			ans += GetZb(i, j) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));
		}
		ans += GetZ1(i, j - 1);
		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	RealScalar SimpleMcCaskillAutoDiff<RealScalar>::GetZb(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);
		if (i == j)
			return RealScalar(0.0);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span)
		{
			return RealScalar(0.0);
		}

		if (Zb[at(i, j)].first)
			return Zb[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);

		// hairpin loop
		ans += exp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, sequence)));
		// internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k)
		{
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l)
			{
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)
					continue;
				ans += GetZb(k, l) * exp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence)));
			}
		}
		// multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k)
		{
			ans += GetZm(i + 1, k - 1) * GetZm1(k, j - 1) * exp(RealScalar(parasor_param::ParDangling(parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]), (j - 1) - 1, (i - 1) + 1, false, sequence))) * exp(RealScalar(parasor_param::ParMultiloopClosing()));
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};

	template <typename RealScalar>
	RealScalar SimpleMcCaskillAutoDiff<RealScalar>::GetZm(const int i, const int j)
	{
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)
			return RealScalar(0.0);
		if (i == j + 1)
			return RealScalar(0.0);
		if (Zm[at(i, j)].first)
			return Zm[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		for (int k = i; k /* + TURN + 1*/ <= j; ++k)
		{
			ans += (1.0 + GetZm(i, k - 1)) * GetZm1(k, j);
		}
		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	RealScalar SimpleMcCaskillAutoDiff<RealScalar>::GetZm1(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);
		if (i == j)
			return RealScalar(0.0);
		if (Zm1[at(i, j)].first)
			return Zm1[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0)
		{
			ans += GetZb(i, j) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence))) * exp(RealScalar(parasor_param::ParMultiloopInternal()));
		}
		ans += GetZm1(i, j - 1);
		return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
	};

	template <typename RealScalar>
	RealScalar SimpleMcCaskillAutoDiff<RealScalar>::GetWb(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span)
		{
			return RealScalar(0.0);
		}

		if (Wb[at(i, j)].first)
			return Wb[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);

		ans += GetZ(1, i - 1) * GetZ(j + 1, n) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));

		//		for (int h = 1; h < i; ++h) {
		for (int h = std::max(1, i - max_span + TURN + 2); h < i; ++h)
		{
			//			for (int l = j + 1; l <= n; ++l) {
			for (int l = j + 1; l <= n && (l - h) <= max_span; ++l)
			{
				const int rtype = parasor_param::GetPairTypeReverse(sequence[h - 1], sequence[l - 1]);
				if (rtype == 0)
					continue;

				const int u = (i - h - 1) + (l - j - 1);
				if (u <= max_loop)
				{
					// internal loop, bulge, stem
					ans += GetWb(h, l) * exp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, sequence)));
				}

				//(h,l)�ŕ�����multiloop
				ans += exp(RealScalar(parasor_param::ParMultiloopClosing())) * exp(RealScalar(parasor_param::ParMultiloopInternal())) * GetWb(h, l) * exp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, sequence))) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence))) * (GetZm(h + 1, i - 1) + GetZm(j + 1, l - 1) + GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1));
			}
		}
		return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	std::pair<std::vector<std::vector<RealScalar>>, std::vector<std::vector<bool>>> SimpleMcCaskillAutoDiff<RealScalar>::run()
	{

		std::cout << "\n[SimpleMcCaskill forward calculation]\n";
		std::cout << "n = " << this->n << std::endl;
		std::cout << "max_span = " << this->max_span << std::endl;
		auto bppm = std::vector<std::vector<RealScalar>>(n + 1, std::vector<RealScalar>(max_span + 1, RealScalar(0.0)));
		auto isbppm_clipped = std::vector<std::vector<bool>>(n + 1, std::vector<bool>(max_span + 1, false));
		// std::cout << "start calculating bppm\n";
		for (int i = 1; i <= n; ++i)
			for (int j = i; j <= n && j - i <= max_span; ++j)
			{
				bppm[i][j - i] = GetZb(i, j) * GetWb(i, j) / GetZ(1, n);

				if (bppm[i][j - i] < 0.0)
				{
					bppm[i][j - i] = 0.0;
					isbppm_clipped[i][j - i] = true; // bppm is clipped --> no gradients
				}
				if (bppm[i][j - i] > 1.0)
				{
					bppm[i][j - i] = 1.0;
					isbppm_clipped[i][j - i] = true; // bppm is clipped --> no gradients
				}
			}

		// std::cout << "end   calculating bppm\n";

		// return std::make_pair(bppm, GetZ(1, n));
		return std::make_pair(bppm, isbppm_clipped);
	}

	template <typename RealScalar>
	std::vector<std::vector<DiffFloat>> SimpleMcCaskillAutoDiff<RealScalar>::diff(
		const std::vector<std::vector<bool>> &isbppm_clipped)
	{

		this->init6DiffVectors(); // clear 6DiffVectors caching

		// std::cout << "\n[SimpleMcCaskill differential calculation]\n";
		auto dbppm_dp = std::vector<std::vector<DiffFloat>>(n + 1, std::vector<DiffFloat>(max_span + 1, DiffFloat(0.0)));
		for (int i = 1; i <= n; ++i)
			for (int j = i; j <= n && j - i <= max_span; ++j)
			{

				if (isbppm_clipped[i][j - i])
					continue; // bppm is clipped --> no gradients

				// (v(x) * u'(x) - u(x) * v'(x)) / [v(x)]^2
				// u(x) = a(x)b(x) (the product of the first two terms)
				// v(x) = c(x) (the term in the denominator)
				// u'(x) = the derivative of u(x) (which will involve the product rule again)
				// v'(x) = the derivative of v(x)
				// bppm[i][j - i] = GetZb(i, j) * GetWb(i, j) / GetZ(1, n);
				auto u = GetZb(i, j) * GetWb(i, j);
				auto dudp = Get_dZb_dp(i, j) * GetWb(i, j) + GetZb(i, j) * Get_dWb_dp(i, j);
				auto v = GetZ(1, n);
				auto dvdp = Get_dZ_dp(1, n);
				dbppm_dp[i][j - i] = (v * dudp - u * dvdp) / (v * v);
				// std::cout << i << ", " << j - i << ", " << u << ", " << dudp << ", " << v << ", " << dvdp << ", " <<  dbppm_dp[i][j - i] << "\n";
			}

		return dbppm_dp;
	}

	template <typename RealScalar>
	DiffFloat SimpleMcCaskillAutoDiff<RealScalar>::Get_dZ_dp(const int i, const int j)
	{
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)
			return DiffFloat(0.0);
		if (i == j + 1)
			return DiffFloat(0.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (dZ_dp[index].first)
			return dZ_dp[index].second;

		DiffFloat ans = DiffFloat(0.0);
		if (i == 1)
		{

			// no base pair
			//  ans += RealScalar(1.0);

			// there is only one outermost base pair (i, *)
			//  ans += GetZ1(i, j);
			ans += Get_dZ1_dp(i, j);

			// there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
			for (int k = i; k <= j - 1; ++k)
			{
				// ans += GetZ(i, k) * GetZ1(k + 1, j);
				ans += Get_dZ_dp(i, k) * GetZ1(k + 1, j) + GetZ(i, k) * Get_dZ1_dp(k + 1, j);
			}
		}
		else
		{
			// i is not paired
			//  ans += GetZ(i + 1, j);
			ans += Get_dZ_dp(i + 1, j);

			//(i,k) form a base pair
			for (int k = i + TURN + 1; k <= j && k - i <= max_span; ++k)
			{
				const int type = parasor_param::GetPairType(sequence[i - 1], sequence[k - 1]);
				if (type == 0)
					continue;
				// ans += GetZb(i, k)
				// 	* exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence)))
				// 	* GetZ(k + 1, j);
				const auto f = exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence)));
				ans += Get_dZb_dp(i, k) * f * GetZ(k + 1, j) + GetZb(i, k) * 0.0 * GetZ(k + 1, j) + GetZb(i, k) * f * Get_dZ_dp(k + 1, j);
			}
		}
		return (dZ_dp[index] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	DiffFloat SimpleMcCaskillAutoDiff<RealScalar>::Get_dZ1_dp(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);
		if (i == j)
			return DiffFloat(0.0);

		if ((j - i) > max_span)
		{
			// return GetZ1(i, i + max_span);
			return Get_dZ1_dp(i, i + max_span);
		}

		if (dZ1_dp[at(i, j)].first)
			return dZ1_dp[at(i, j)].second;

		DiffFloat ans = DiffFloat(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0)
		{
			const auto f = exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));
			// ans += GetZb(i, j) * f;
			ans += Get_dZb_dp(i, j) * f + GetZb(i, j) * 0.0;
		}
		// ans += GetZ1(i, j - 1);
		ans += Get_dZ1_dp(i, j - 1);

		return (dZ1_dp[at(i, j)] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	DiffFloat SimpleMcCaskillAutoDiff<RealScalar>::Get_dZb_dp(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);
		if (i == j)
			return DiffFloat(0.0);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span)
		{
			return DiffFloat(0.0);
		}

		if (dZb_dp[at(i, j)].first){
			return dZb_dp[at(i, j)].second;
		}
		
		DiffFloat ans = DiffFloat(0.0);

		// hairpin loop
		//  ans += exp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, sequence)));
		ans += DiffFloat(0.0);

		// internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k)
		{
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l)
			{
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)
					continue;
				const auto f = exp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence)));
				const auto dfdp = f * parasor_param::dParLoopEnergy_dp(i - 1, j - 1, k - 1, l - 1, sequence, pidx1, pidx2);
				// ans += GetZb(k, l) * f;
				ans += Get_dZb_dp(k, l) * f + GetZb(k, l) * dfdp;
			}
		}

		// multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k)
		{
			const auto f1 = exp(RealScalar(parasor_param::ParDangling(parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]), (j - 1) - 1, (i - 1) + 1, false, sequence)));
			const auto f2 = exp(RealScalar(parasor_param::ParMultiloopClosing()));
			// ans += GetZm(i + 1, k - 1) * GetZm1(k, j - 1) * f1 * f2;
			ans += Get_dZm_dp(i + 1, k - 1) * GetZm1(k, j - 1) * (f1 * f2) + GetZm(i + 1, k - 1) * Get_dZm1_dp(k, j - 1) * (f1 * f2) + GetZm(i + 1, k - 1) * GetZm1(k, j - 1) * (0.0);
		}

		return (dZb_dp[at(i, j)] = std::make_pair(true, ans)).second;
	};

	template <typename RealScalar>
	DiffFloat SimpleMcCaskillAutoDiff<RealScalar>::Get_dZm_dp(const int i, const int j)
	{
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)
			return DiffFloat(0.0);
		if (i == j + 1)
			return DiffFloat(0.0);
		if (dZm_dp[at(i, j)].first)
			return dZm_dp[at(i, j)].second;

		DiffFloat ans = DiffFloat(0.0);
		for (int k = i; k /* + TURN + 1*/ <= j; ++k)
		{
			// ans += (1.0 + GetZm(i, k - 1)) * GetZm1(k, j);
			ans += (0.0 + Get_dZm_dp(i, k - 1)) * GetZm1(k, j) + (1.0 + GetZm(i, k - 1)) * Get_dZm1_dp(k, j);
		}
		return (dZm_dp[at(i, j)] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	DiffFloat SimpleMcCaskillAutoDiff<RealScalar>::Get_dZm1_dp(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);
		if (i == j)
			return DiffFloat(0.0);
		if (dZm1_dp[at(i, j)].first)
			return dZm1_dp[at(i, j)].second;

		DiffFloat ans = DiffFloat(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0)
		{
			const auto f1 = exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)));
			const auto f2 = exp(RealScalar(parasor_param::ParMultiloopInternal()));
			// ans += GetZb(i, j) * f1 * f2;
			ans += Get_dZb_dp(i, j) * (f1 * f2) + GetZb(i, j) * (0.0);
		}
		// ans += GetZm1(i, j - 1);
		ans += Get_dZm1_dp(i, j - 1);
		return (dZm1_dp[at(i, j)] = std::make_pair(true, ans)).second;
	};

	template <typename RealScalar>
	DiffFloat SimpleMcCaskillAutoDiff<RealScalar>::Get_dWb_dp(const int i, const int j)
	{
		assert(1 <= i && i <= j && j <= n);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span)
		{
			return DiffFloat(0.0);
		}

		if (dWb_dp[at(i, j)].first)
			return dWb_dp[at(i, j)].second;

		DiffFloat ans = DiffFloat(0.0);

		const auto f1 = exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));
		// ans += GetZ(1, i - 1) * GetZ(j + 1, n) * f1;
		ans += Get_dZ_dp(1, i - 1) * GetZ(j + 1, n) * f1 + GetZ(1, i - 1) * Get_dZ_dp(j + 1, n) * f1 + GetZ(1, i - 1) * GetZ(j + 1, n) * 0.0;

		//		for (int h = 1; h < i; ++h) {
		for (int h = std::max(1, i - max_span + TURN + 2); h < i; ++h)
		{
			//			for (int l = j + 1; l <= n; ++l) {
			for (int l = j + 1; l <= n && (l - h) <= max_span; ++l)
			{
				const int rtype = parasor_param::GetPairTypeReverse(sequence[h - 1], sequence[l - 1]);
				if (rtype == 0)
					continue;

				const int u = (i - h - 1) + (l - j - 1);
				if (u <= max_loop)
				{
					// internal loop, bulge, stem
					const auto f = exp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, sequence)));
					const auto dfdp = f * parasor_param::dParLoopEnergy_dp(h - 1, l - 1, i - 1, j - 1, sequence, pidx1, pidx2);
					// ans += GetWb(h, l) * f;
					ans += Get_dWb_dp(h, l) * f + GetWb(h, l) * dfdp;
				}

				//(h,l)�ŕ�����multiloop
				const auto f2 = exp(RealScalar(parasor_param::ParMultiloopClosing())) * exp(RealScalar(parasor_param::ParMultiloopInternal())) * exp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, sequence))) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)));
				// ans += f2 * GetWb(h, l) * ( GetZm(h + 1, i - 1) + GetZm(j + 1, l - 1) + GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1) );
				const auto f3 = GetZm(h + 1, i - 1) + GetZm(j + 1, l - 1) + GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1);
				const auto df3dp = Get_dZm_dp(h + 1, i - 1) + Get_dZm_dp(j + 1, l - 1) + Get_dZm_dp(h + 1, i - 1) * GetZm(j + 1, l - 1) + GetZm(h + 1, i - 1) * Get_dZm_dp(j + 1, l - 1);
				ans += 0.0 * GetWb(h, l) * f3 + f2 * Get_dWb_dp(h, l) * f3 + f2 * GetWb(h, l) * df3dp;
			}
		}
		return (dWb_dp[at(i, j)] = std::make_pair(true, ans)).second;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::setSequence(const std::string sequence, const int max_span)
	{
		this->sequence = sequence;
		this->max_span = max_span;

		this->n = int(sequence.size());
		this->init6Vectors();
		std::cout << "set sequence to \n " << this->sequence << "\n";
		std::cout << "set max_span to : " << this->max_span << "\n";

		// this->init6DiffVectors();  // perform this initilaization in the train method
		return;
	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::setPidx(const int pidx1, const int pidx2)
	{
		this->pidx1 = pidx1;
		this->pidx2 = pidx2;
		return;
	}

	template <typename RealScalar>
	std::vector<RealScalar> SimpleMcCaskillAutoDiff<RealScalar>::calc_bppmsum(
		const std::vector<std::vector<RealScalar>> &bppm,
		const std::string seq,
		const int max_span)
	{
		// std::cout << bppm.size() << "\n";
		// std::cout << bppm[0].size()<< "\n";
		std::vector<RealScalar> bppm_sum(seq.size(), RealScalar(0.0));
		const int max_span_constraint = max_span;
		for (int i = 0; i < int(seq.size()); ++i)
		{
			int upper = std::max(0, i - max_span_constraint);
			for (int j = i; j < upper; ++j)
			{

				bppm_sum[i] += bppm[j + 1][i - j];
			}
			for (int j = 0; j < max_span_constraint + 1; ++j)
			{
				bppm_sum[i] += bppm[i + 1][j];
			}
		}

		return bppm_sum;
	}

	template <typename RealScalar>
	std::vector<std::pair<std::string, typename SimpleMcCaskillAutoDiff<RealScalar>::BppmType>> SimpleMcCaskillAutoDiff<RealScalar>::load_sequence_and_bppm(
		const std::string datafile)
	{
		
		std::vector<std::pair<std::string, BppmType>> seqbppms;
		BppmType bppm;

		std::ifstream bppmfile(datafile);
		if (!bppmfile.is_open()) {
			std::cerr << "Unable to open file :\n  " << datafile << std::endl;
			return seqbppms;
		}

		

		std::string line;
		std::string seq("");
		while (std::getline(bppmfile, line)) {
			// skip empty row

			if (line.empty()) continue;

			
			// read sequence
			if (std::isalpha(line[0])) {
				if ((seq.size()>0) && (bppm.size()>0)) {
					bppm.insert(bppm.begin(), std::vector<RealScalar>(bppm.at(0).size() + 1, RealScalar(0.0)));
					seqbppms.push_back(std::make_pair(seq, bppm));
					
					// std::cout << "added item to seqbppms : " << seqbppms.size() << std::endl;
					std::cout << "loaded sequences : " << seq << std::endl;
					//debug
					std::cout << "loaded bppm size: " << bppm.size() << "&" << bppm[0].size() << std::endl;
					for (size_t i = 0; i < bppm.size(); i++){
						std::cout << "bppm[" << i << "]:size" << bppm[i].size() << std::endl;
						for (size_t j = 0; j < bppm[i].size(); j++){
							std::cout << bppm[i][j] << " ";
						}
						std::cout << std::endl;
					}
					//debug end
				}
				bppm.clear();
				seq = trim(line);
				continue;		
			}

			if (std::isdigit(line[0])) {
				std::stringstream ss(line);
				std::string value;
				std::vector<RealScalar> values;
				values.push_back(RealScalar(0.0));
				while (std::getline(ss, value, ',')) {
					values.push_back(std::stod(value));
				}
				bppm.push_back(values);
				//debug 50配列渡すと、valueには51個の数値が入る（0.0をpush_backしてるから）
				#if 0
				std::cout << "values total_size:" << values.size() << std::endl;
				for (size_t i = 0; i < values.size(); ++i) {
        			std::cout << std::fixed << std::setprecision(6) << values[i] << " ";
    			}
    			std::cout << std::endl;
				std::cout << "size bppm:" << bppm.size() << std::endl;
				#endif
				//debug end
			}
		}
		if ((seq.size()>0) && (bppm.size()>0)) {
			//この処理が謎。サイズ+1の0.0の配列をインサートしている
			//最後の配列のbppmしか残していない？
			bppm.insert(bppm.begin(), std::vector<RealScalar>(bppm.at(0).size() + 1, RealScalar(0.0)));
			seqbppms.push_back(std::make_pair(seq, bppm));
			std::cout << "loaded sequences : " << seq << std::endl;
			//bppmは50塩基が入ると51x52の配列ができる（謎）
			//debug
			std::cout << "loaded bppm size: " << bppm.size() << "&" << bppm[0].size() << std::endl;
			for (size_t i = 0; i < bppm.size(); i++){
				std::cout << "bppm[" << i << "]:size" << bppm[i].size() << std::endl;
				for (size_t j = 0; j < bppm[i].size(); j++){
					std::cout << bppm[i][j] << " ";
				}
				std::cout << std::endl;
			}
			//debug end
			bppm.clear();
			seq = "";
		}

		bppmfile.close();

		return seqbppms;
	}

	template <typename RealScalar>
	std::vector<DiffFloat> SimpleMcCaskillAutoDiff<RealScalar>::calc_dbppmsum(
		const std::vector<std::vector<DiffFloat>> &dbppm_dp,
		const std::string seq,
		const int max_span)
	{
		// std::cout << bppm.size() << "\n";
		// std::cout << bppm[0].size()<< "\n";
		std::vector<DiffFloat> dbppm_sum(seq.size(), DiffFloat(0.0));
		const int max_span_constraint = max_span;
		for (int i = 0; i < int(seq.size()); ++i)
		{
			int upper = std::max(0, i - max_span_constraint);
			for (int j = i; j < upper; ++j)
			{

				dbppm_sum[i] += dbppm_dp[j + 1][i - j];
			}
			for (int j = 0; j < max_span_constraint + 1; ++j)
			{
				dbppm_sum[i] += dbppm_dp[i + 1][j];
			}
		}

		return dbppm_sum;
	}


	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::calc(
				const std::string datafile,
				const std::string outputfile
	)
	{
		std::vector<std::string> sequences;

		/* load update mask */
		this->load_update_mask("update_mask.txt");
		
		std::cout << "loading sequences from :\n  " << datafile << std::endl;
		/* read m6a dataset */
		load_m6a_dataset(datafile, sequences);

		/* to log bppm */
		std::string bppmfile;
		if(outputfile.empty()){
			bppmfile = replace_string(datafile, ".fa", ".bppm");
		}else{
			bppmfile = outputfile;
		}

		/* calculate bppmsum as answer */
		for (int j=0; j < sequences.size(); ++j) {
			std::cout << "\ntraining data No." << j + 1 << "\n";
			auto& seq = sequences[j];
			this->setSequence(seq, seq.size());
			/*bppmはファイルからロードされていればそれを使う、ロードされていなければ計算する */
			const std::vector<std::vector<RealScalar>> bppm = this->run().first;
			this->write_bppm(bppmfile, bppm, seq);
		}


	}

	template <typename RealScalar>
	void SimpleMcCaskillAutoDiff<RealScalar>::train(
		const int nepoch,
		const double learning_rate,
		int nsample,
		const int sequence_length,
		int max_span,
		const std::string logfile,
		const double logstack_init,
		const double weight_central,
		const std::string datafile)
	{

		/* load update mask */
		this->load_update_mask("update_mask.txt");

		std::ofstream outfile(logfile, std::ios::app);
		if (outfile.is_open())
		{
			outfile << "[training]\n";
		}

		outfile.close();

		/*
		 * step 1. make training data
		 */
		std::cout << "[make training data]\n";

		std::vector<std::string> sequences;
		std::vector<std::vector<RealScalar>> bppmsums;
		
		parasor_param::print_logstack(logfile, this->update_mask);


		/* makeing random sequence */
		bool is_bppm_loaded = false;
		std::vector<std::pair<std::string, BppmType>> seqbppms;
		if (datafile.empty()) {
			std::cout << "generating random sequences\n";
			for (int j = 0; j < nsample; ++j)
			{
				// const std::string lfile = (j == 0) ? logfile : "";
				// parasor_param::print_logstack(lfile, this->update_mask);

				


				std::cout << "training data No." << j + 1 << "\n";

				auto seq = MakeRandomRnaSeq(sequence_length, sequence_length + j); // (length, seed)
				sequences.push_back(seq);

				/* moved to indivisual loop */
				// this->setSequence(seq, max_span);
				// const std::vector<std::vector<RealScalar>> bppm = this->run().first;
				// const std::vector<RealScalar> bppmsum = this->calc_bppmsum(bppm, seq, max_span);
				// bppmsums.push_back(bppmsum);
				// std::cout << "bppmsum :\n " << bppmsum[0] << ", " << bppmsum[1] << ", " << bppmsum[2] << " ..." << "\n";
			}
		}
		/* loading sequence from dataset */
		else {
			if (endsWith(datafile, ".fa")) {
				std::cout << "loading sequences from :\n  " << datafile << std::endl;
				/* read m6a dataset */
				load_m6a_dataset(datafile, sequences);
			}
			else if (datafile, ".bppm") {
				std::cout << "loading sequences and bppms from :\n  " << datafile << std::endl;
				seqbppms = this->load_sequence_and_bppm(datafile);
				if (seqbppms.size()==0) {
					std::cerr << "failed to read sequences and bppm data from :\n  " << datafile << std::endl;
					return;
				}

				max_span = seqbppms.at(0).first.size();
				std::cout << "max_span is set to : " << max_span << std::endl;
				for (const auto& elem : seqbppms) {
					auto& seq = elem.first;
					auto& bppm = elem.second;
					sequences.push_back(seq);
				}
				is_bppm_loaded = true;
			}
			else {
				std::cerr << "unknown extension of file :\n " << datafile << std::endl;
				return;
			}

			/* read m6a dataset */
			nsample = sequences.size();
			if (nsample < 1) {
				std::cerr << "nsample = " << nsample << std::endl;
				return;
			}
		}		
		std::cout << "number of training data : " << sequences.size() << std::endl;
		std::cout << "length of each sequence :\n";
		for (auto &seq : sequences)
		{
			std::cout << seq.size() << ", ";
		}
		std::cout << "\n";

		/* to log bppm */
		std::string bppmfile = replace_string(logfile, ".txt", ".bppm");
		
		/* calculate bppmsum as answer */
		for (int j=0; j < sequences.size(); ++j) {
			std::cout << "\ntraining data No." << j + 1 << "\n";
			auto& seq = sequences[j];
			this->setSequence(seq, max_span);
			/*bppmはファイルからロードされていればそれを使う、ロードされていなければ計算する */
			const std::vector<std::vector<RealScalar>> bppm = is_bppm_loaded? seqbppms.at(j).second : this->run().first;
			this->write_bppm(bppmfile, bppm, seq);
			const std::vector<RealScalar> bppmsum = this->calc_bppmsum(bppm, seq, max_span);
			bppmsums.push_back(bppmsum);
			std::cout << "bppmsum :\n " << bppmsum[0] << ", " << bppmsum[1] << ", " << bppmsum[2] << " ..." << "\n";
		}


		/*
		 * step 2. setup training, start training epoch
		 */

		std::cout << "[training]\n";
		// set logstack parameter to 0 as training initialization
		parasor_param::reset_logstack(logstack_init, weight_central, this->update_mask);
		parasor_param::update_parstack();
		for (int epoch = 1; epoch <= nepoch; ++epoch)
		{
			std::cout << "[epoch = " << epoch << "]\n";
			std::ofstream outfile(logfile, std::ios::app);
			if (outfile.is_open())
			{
				outfile << "\n[epoch = " << epoch << "]\n";
			}
			outfile.close();

			for (int i = 0; i < nsample; ++i)
			{	// loop thorugh sequences
				/*
				 * step 3. forward calculation for bmmp
				 */
				std::cout << "  [start training sequence]\n";

				const auto start = std::chrono::system_clock::now();

				auto &seq = sequences[i];
				this->setSequence(seq, max_span);
				auto ret = this->run();
				auto bppm = ret.first;
				auto isbppm_clipped = ret.second;
				auto bppmsum = this->calc_bppmsum(bppm, seq, max_span);
				auto &bppmsum_answer = bppmsums[i];

				/*
				 * step 4. compute gradients with respective to each single parameter in logstack[13][13], and store it
				 */
				// loop through logstack 2d index, pidx1, pidx2
				for (int k1 = 0; k1 < this->pidx1Size; ++k1)
				{
					for (int k2 = 0; k2 < this->pidx2Size; ++k2)
					{
						if (!this->update_mask[k1][k2]) continue;

						// std::cout << "(k1,k2) = (" << k1 << "," << k2 << ")\n";
						this->setPidx(k1 + this->pidx1_lower, k2 + this->pidx2_lower);
						auto dbppm_dp = this->diff(isbppm_clipped);
						// print2DVector<DiffFloat>(dbppm_dp);

						auto dbppmsum_dp = this->calc_dbppmsum(dbppm_dp, seq, max_span);
						DiffFloat dp = 0.0;
						// const double fac = (((k1==2)&&(k2==2)) || ((k1==3)&&(k2==3)))? 1.0 : weight_central;
						const double fac = 1.0;
						for (int k = 0; k < bppmsum.size(); ++k)
						{
							auto error = bppmsum[k] - bppmsum_answer[k];
							//変更 signの符号が逆
							//double sign = (error > 0) - (error < 0);
							double sign = (error < 0) - (error > 0);
							// std::cout << learning_rate << ", " << error << ","  << dbppmsum_dp[k] << std::endl;
							// dp += - learning_rate * error * dbppmsum_dp[k] * fac;   // L2Loss
							dp += -learning_rate * sign * dbppmsum_dp[k] * fac; // L1Loss
							
							//debug
							#if 0
							if ( (k1 == 0) & (k2 == 10)){
								std::cout << "dp:" << dp << ",error : " << error << ", sign: " << sign << ", learning_rate: " << learning_rate << ", dbppmsum_dp: " << dbppmsum_dp[k] << ",k:" << k << ",fac:" << fac << std::endl;
							}
							#endif
						}
						this->add_dlogstack(dp, k1, k2);
						//debug
						#if 0
						if ( (k1 == 0) & (k2 == 10)) {
							std::cout << "last_dp:" << dp << ",k1:" << k1 << ",k2:" << k2 << ",:stack:" << this->dlogstack[k1][k2] << std::endl;
						}
						#endif
					}
				}

				// this->print_dlogstack();
				/*
				* step 4.5 make symmetry dlogstack 
				*/
				for (int k1 = 0; k1 < this->pidx1Size; ++k1) for (int k2 = k1+1; k2 < this->pidx2Size; ++k2) {
					this->dlogstack[k1][k2] += this->dlogstack[k2][k1];
					this->dlogstack[k1][k2] *= 0.5;
					this->dlogstack[k2][k1] = this->dlogstack[k1][k2];
				}

				// this->print_dlogstack();	


				// this->print_dlogstack();
				/*
				 * step 5. updates logstack[13][13]
				 */
				for (int k1 = 0; k1 < this->pidx1Size; ++k1)
				{
					for (int k2 = 0; k2 < this->pidx2Size; ++k2)
					{
						/*
						* although this->dlogstack[k1][k2] is 0 if this->update_mask[k1][k2]==false
						* for the sake of safety we skip the update of 
						* logstack[k1 + this->pidx1_lower][k2 + this->pidx2_lower]
						*/
						if (!this->update_mask[k1][k2]) continue;
						parasor_param::update_logstack(this->dlogstack[k1][k2], k1 + this->pidx1_lower, k2 + this->pidx2_lower);
					}
				}
				parasor_param::update_parstack();
				this->reset_dlogstack(); // reset dlogstack every sequence --> bacthsize=1
				
				const std::string lfile = (i == nsample - 1) ? logfile : "";
				parasor_param::print_logstack(lfile, this->update_mask);

				const auto end = std::chrono::system_clock::now();
				const int time = int(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
				std::cout << "( forward + backward ) execution time : " << time << " [msec]\n";
				std::cout << "[epoch = " << epoch << "]\n";
			}
		}
	}

}
