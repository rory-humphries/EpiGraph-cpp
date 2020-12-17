//
// Created by roryh on 22/06/2020.
//
#include <EpiGraph/Random/Distributions.hpp>

#include <random>

namespace EpiGraph {
	std::mt19937 &global_engine() {
		static std::random_device rd{};
		static std::mt19937 gen(rd());
		//static std::default_random_engine e{};
		return gen;
	}

	auto ProbDist_from_csv(std::string csv) -> ProbDist {
		std::ifstream infile;
		infile.open(csv);
		if (infile.fail()) {
			throw std::runtime_error("Failed to open file " + csv);
		}
		std::string line;

		std::vector<std::pair<double, double>> val_prob_pairs;

		double val;
		double prob;

		getline(infile, line, '\n'); // header

		while (getline(infile, line, '\n')) {
			std::istringstream ss(line);
			std::string token;

			std::getline(ss, token, ',');

			val = stod(token);

			std::getline(ss, token, ',');

			prob = stod(token);

			val_prob_pairs.emplace_back(val, prob);
		}
		infile.close();
		return ProbDist(val_prob_pairs.begin(), val_prob_pairs.end());

	}

	ProbDist &ProbDist::operator=(const ProbDist &rhs) {
		index_dist = rhs.index_dist;
		vals = rhs.vals;
		probs = index_dist.probabilities();

		return *this;
	}


	auto ProbDist::get_prob(double val) -> double {
		auto it = std::lower_bound(vals.begin(), vals.end(), val);
		if (it == vals.end())
			it--;
		int index = std::distance(vals.begin(), it);

		return index_dist.probabilities().at(index);
	}
}