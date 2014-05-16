
#include <iostream>
#include <ctime>
#include <cmath>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"

int main()
{

    gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    int seed = time(NULL);
    gsl_rng_set(rng, seed);

    double var1 = 25.0;
    double var2 = 1000.0;

    double total_var = var1 + var2;

    double a = sqrt(3.0*total_var);

    std::cerr << a << std::endl;


    double x1_mean = 0.0;
    double x1_std = 0.0;

    double x2_mean = 0.0;
    double x2_std = 0.0;

    double s_mean = 0.0;
    double s_std = 0.0;

    int N = 1000;

    int M = 10000;

    for(int i = 0; i < N; i++)
    {
        double x1 = gsl_ran_gaussian(rng, sqrt(var1));

        double x2 = 0.0;

        for(int j = 0; j < M; j++)
        {
            double s = gsl_ran_gaussian(rng, sqrt(var2)) + x1;

            double cdf = gsl_cdf_gaussian_P(s, sqrt(total_var));

            //s = gsl_cdf_flat_Pinv(cdf, -a, a);
            //s = gsl_cdf_gaussian_Pinv(cdf, sqrt(total_sigma));

            double lognormal_sigma = sqrt(log((1.0+sqrt(1.0+4*total_var))/2.0));

            s = gsl_cdf_lognormal_Pinv(cdf, 0.0, lognormal_sigma);

            x2 = x2 + s;

            s_mean += s;
            s_std += s * s;
        }



        x2 = x2 / 10000;

        std::cout << x1 << " " << x2 << std::endl;

        x1_mean += x1;
        x2_mean += x2;

        x1_std += x1 * x1;
        x2_std += x2 * x2;
    }

    x1_mean /= N;
    x2_mean /= N;

    x1_std /= N;
    x2_std /= N;

    s_mean /= (N*M);
    s_std /= (N*M);

    x1_std -= x1_mean * x1_mean;
    x2_std -= x2_mean * x2_mean;

    std::cerr << "x1_mean: " << x1_mean << " x1_std: " << x1_std << std::endl;
    std::cerr << "x2_mean: " << x2_mean << " x2_std: " << x2_std << std::endl;
    std::cerr << "s_mean: " << s_mean << " s_std: " << s_std << std::endl;




    //std::cout << x1 << " " << net <<  std::endl;

    return 0;
}
