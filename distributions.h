//
// Created by roryh on 05/04/2020.
//

#ifndef EPIGRAPH_DISTRIBUTIONS_H
#define EPIGRAPH_DISTRIBUTIONS_H

#include <cmath>
#include <random>
#include <fstream>
#include <sstream>

// not thread safe, careful with openMP
std::mt19937 &global_engine() {
    static std::random_device rd{};
    static std::mt19937 gen(rd());
    //static std::default_random_engine e{};
    return gen;
}

class power_law{
    /*
     * Models a power law distribution.
     */

public:
    std::uniform_real_distribution<> uni_dist;
    double m_alpha;
    double m_xmin;
    double m_xmax;

    power_law(double alpha, double xmin): m_alpha(alpha), m_xmin(xmin),  uni_dist(0.0, 1.0) {}

    static double pdf(double x, double alpha, double x_min)
    {
        /*
         * The pdf for a power law distribution
         */
        return ((alpha - 1) / x_min) * pow((x / x_min) , -alpha);
    }
    static double non_normal_pdf(double x, double alpha, double x_min)
    {
        /*
         * The non-normalised pdf (does not integrate to 1) for a power law distribution
         */
        return pow(x, -alpha) ;
    }
    static double cdf(double u, double alpha, double x_min)
    {
        /*
         * The cdf for a power law distribution
         */
        return pow(u / x_min, -alpha + 1);
    }
    static double cdf_inv(double u, double alpha, double x_min)
    {
        /*
        * The inverse cdf for a power law distribution
        */
        return pow(x_min * u, -1 / (alpha - 1));
    }
    template<typename RandomGenerator>
    double operator()( RandomGenerator& gen) {
        double r = uni_dist(gen);
        return cdf_inv(r, m_alpha, m_xmin);
    }
};

struct prob_dist
{
    std::vector<double> cumulative_dist;
    std::vector<double> values;
    std::uniform_real_distribution<double> uni01_distribution;

    prob_dist()= default;
    explicit prob_dist(std::string csv)
    {

        std::uniform_real_distribution<double> uni01_distribution_tmp(0.0,1.0);
        uni01_distribution = uni01_distribution_tmp;

        double sum = 0;

        std::ifstream myfile;
        std::string line;
        myfile.open(csv);

        if (myfile.is_open())
        {
            getline(myfile, line);

            while (getline(myfile, line))
            {
                std::istringstream s(line);
                std::string field;

                while (getline(s, field, ','))
                {
                    sum += stod(field);
                    cumulative_dist.push_back(sum);

                    getline(s, field, ',');
                    values.push_back(stod(field));
                }
            }
        }
    }

    template<typename RandomGenerator>
    double operator()( RandomGenerator& gen) {

        double r = uni01_distribution(gen);
        for (int i = 0; i < cumulative_dist.size(); i++)
        {
            if (r < cumulative_dist[i])
            {
                return values[i];
            }
        }
        return r;
    }
};

#endif //EPIGRAPH_DISTRIBUTIONS_H
