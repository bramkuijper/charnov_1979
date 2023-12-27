#include <vector>
#include <iostream>
#include <fstream>
#include <random>

// random number generators
std::random_device rd;
unsigned seed = rd();
std::mt19937 rng_r(seed);

// parameters
unsigned max_generation = 10000;

// mutation rates
double mu_r = 0.02;
double sdmu_r = 0.02;

// for stats purposes
double meanw_gain_curvem = 0.0;
double meanw_gain_curvef = 0.0;
double meanwf = 0.0;
double meanwm = 0.0;

// uniform dist between 0 and 1
std::uniform_real_distribution<double> uniform{0.0,1.0};
// normal dist between 0 and 1
std::normal_distribution<double> mutational_effects{sdmu_r};

unsigned popsize = 5000;

double n = 0.8;
double r_init = 0.1;



std::vector <double> r(popsize,r_init);

std::string file_name = "data.csv";

std::ofstream data_file(file_name);

void mutate() 
{
    for (unsigned idx = 0; idx < popsize; ++idx)
    {
        if (uniform(rng_r) < mu_r)
        {
            r[idx] = std::clamp(r[idx] + mutational_effects(rng_r), 0.0, 1.0);
        }
    }
}

void reproduce()
{
    //allocate vectors of gain curves 
    std::vector<double> wf(popsize,0.0);
    std::vector<double> wm(popsize,0.0);

    meanw_gain_curvem = 0.0;
    meanw_gain_curvef = 0.0;

    // calculate gain curves and stats
    for (unsigned idx = 0; idx < popsize; ++idx)
    {
        wf[idx] = 1.0 - r[idx];

        meanw_gain_curvef += wf[idx];

        wm[idx] = std::pow(r[idx], n);

        meanw_gain_curvem += wm[idx];
    }

    meanw_gain_curvem /= popsize;
    meanw_gain_curvef /= popsize;

    std::vector r_new_generation(popsize, 0.0);

    // make sampling distribution of maternal gametes
    // (i.e., sampling according to fitness weightings
    std::discrete_distribution<unsigned> maternal_allele_sampler(wf.begin(), wf.end());
    std::discrete_distribution<unsigned> paternal_allele_sampler(wm.begin(), wm.end());

    meanwf = 0.0;
    meanwm = 0.0;

    for (unsigned idx = 0; idx < popsize; ++idx)
    {
        if (uniform(rng_r) < 0.5)
        {
            // meiosis determines we need a paternal allele
            r_new_generation[idx] = r[paternal_allele_sampler(rng_r)];

            ++meanwm;
        }
        else
        {
            // meiosis determines we need a paternal allele
            r_new_generation[idx] = r[maternal_allele_sampler(rng_r)];
            ++meanwf;
        }
    }

    r = r_new_generation;
}

void write_data(unsigned const time_step)
{
    double meanr = 0.0;
    double varr = 0.0;

    for (unsigned idx = 0; idx < popsize; ++idx)
    {
        meanr += r[idx];
        varr += r[idx] * r[idx];
    }

    meanr /= r.size();
    varr = varr/r.size() - meanr * meanr;

    data_file << time_step << "," << meanr << "," << varr << "," << meanw_gain_curvef << "," << meanw_gain_curvem << "," << meanwm << "," << meanwf << std::endl;
}

void write_header()
{
    data_file << "time,meanr,varr,meanw_gain_curvef,meanw_gain_curvem,meanwf,meanwm" << std::endl;
}

int main()
{
    write_header();

    for (unsigned generation = 0; generation < max_generation; ++generation)
    {
        write_data(generation);
        mutate();
        reproduce();

    }

}
