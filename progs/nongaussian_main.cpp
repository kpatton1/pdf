
#include <math.h>
#include "Field.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_rng.h"

#include "CDF.h"

#include <ctime>


#include <iostream>
#include <fstream>

//#define BREAKPOINT_K 0.1

#define GRID_N_MACRO 128
#define GRID_N_MICRO 16
//#define GRID_L_MACRO (2.0 * M_PI * GRID_N / BREAKPOINT_K)
#define GRID_L_MACRO 119.0
#define GRID_L_MICRO (GRID_L_MACRO/GRID_N_MACRO)

#define ALPHA -1.0

#define PDF_SIZE 99
#define PS_SIZE 100000

double k_array[PS_SIZE];
double pk_array[PS_SIZE];

double l_array_macro[PS_SIZE+1];
double pl_array_macro[PS_SIZE+1];

double l_array_micro[PS_SIZE+1];
double pl_array_micro[PS_SIZE+1];

gsl_interp* pl_interp_macro;
gsl_interp_accel* pl_interp_macro_accel;

gsl_interp* pl_interp_micro;
gsl_interp_accel* pl_interp_micro_accel;

CDF *cdf;

double macro_pl(double l)
{
    if( l <= 0.0)
    {
        return 0.0;
    }

    if(l < l_array_macro[0])
    {
        std::cerr << "macro_pl evaluation low! adjusting l: " << l << " -> " << l_array_macro[0] << std::endl;
        l = l_array_macro[0];
    }
    if(l > l_array_macro[PS_SIZE])
    {
        std::cerr << "macro_pl evaluation high! adjusting l: " << l << " -> " << l_array_macro[PS_SIZE] << std::endl;
        l = l_array_macro[PS_SIZE];
    }

    double pl = gsl_interp_eval(pl_interp_macro, l_array_macro, pl_array_macro, l, pl_interp_macro_accel);

    return pl;
}

double micro_pl(double l)
{
    if( l <= 0.0)
    {
        return 0.0;
    }

    if(l < l_array_micro[0])
    {
        std::cerr << "micro_pl evaluation low! adjusting l: " << l << " -> " << l_array_micro[0] << std::endl;
        l = l_array_micro[0];
    }
    if(l > l_array_micro[PS_SIZE])
    {
        std::cerr << "micro_pl evaluation high! adjusting l: " << l << " -> " << l_array_micro[PS_SIZE] << std::endl;
        l = l_array_micro[PS_SIZE];
    }

    double pl = gsl_interp_eval(pl_interp_micro, l_array_micro, pl_array_micro, l, pl_interp_micro_accel);


    return pl;
}

void init()
{


    double lnx_array[PDF_SIZE];
    double pdf_array[PDF_SIZE];

    std::ifstream ps_file("ps.dat", std::ifstream::in);
    std::ifstream pdf_file("pdf.dat", std::ifstream::in);

    l_array_macro[0] = 0.0;
    pl_array_macro[0] = 0.0;

    l_array_micro[0] = 0.0;
    pl_array_micro[0] = 0.0;

    for(int i = 0; i < PS_SIZE; i++)
    {
        if(!ps_file.good())
        {
            std::cerr << "Error reading ps.dat!" << std::endl;
            break;
        }

        ps_file >> k_array[i] >> pk_array[i];

        l_array_macro[i+1] = GRID_L_MACRO * k_array[i] / (2.0 * M_PI);
        l_array_micro[i+1] = GRID_L_MICRO * k_array[i] / (2.0 * M_PI);

        pl_array_macro[i+1] = pk_array[i] / (2.0 * M_PI  * l_array_macro[i+1] * l_array_macro[i+1]);
        pl_array_micro[i+1] = pk_array[i] / (2.0 * M_PI  * l_array_micro[i+1] * l_array_micro[i+1]);
    }

    for(int i = 0; i < PDF_SIZE; i++)
    {
        if(!pdf_file.good())
        {
            std::cerr << "Error reading ps.dat!" << std::endl;
            break;
        }

        pdf_file >> lnx_array[i] >> pdf_array[i];
    }

    std::cerr << "Loaded files!" << std::endl;

    pl_interp_macro = gsl_interp_alloc(gsl_interp_cspline, PS_SIZE+1);
    pl_interp_macro_accel = gsl_interp_accel_alloc();
    gsl_interp_init(pl_interp_macro, l_array_macro, pl_array_macro, PS_SIZE+1);

    pl_interp_micro = gsl_interp_alloc(gsl_interp_cspline, PS_SIZE+1);
    pl_interp_micro_accel = gsl_interp_accel_alloc();
    gsl_interp_init(pl_interp_micro, l_array_micro, pl_array_micro, PS_SIZE+1);

    std::cerr << "GRID_N_MACRO: " << GRID_N_MACRO << "  GRID_N_MICRO: " << GRID_N_MICRO << "  GRID_L_MACRO: " << GRID_L_MACRO << "  GRID_L_MICRO: " << GRID_L_MICRO << std::endl;
    std::cerr << "macro grid: kmin = " << (2.0*M_PI*1.0 / GRID_L_MACRO) << "  kmax = " << (2.0*M_PI*sqrt(2.0)*GRID_N_MACRO / 2.0 / GRID_L_MACRO) << std::endl;
    std::cerr << "micro grid: kmin = " << (2.0*M_PI*1.0 / GRID_L_MICRO) << "  kmax = " << (2.0*M_PI*sqrt(2.0)*GRID_N_MICRO / 2.0 / GRID_L_MICRO) << std::endl;

    std::cerr << "Set up ps interp!" << std::endl;

    cdf = new CDF_PDF_INTERP(lnx_array, pdf_array, PDF_SIZE);

    std::cerr << "Set up CDF!" << std::endl;
}



/*void print(double l, double pl)
{
    double k = 2.0 * M_PI * l / GRID_L;
    double pk = (2.0 * M_PI * l * l) * pl;
    
    std::cout << k << " " << pk << std::endl;
}*/

int main(int argc, char* argv[])
{
    init();

    Field2d macro_field(GRID_N_MACRO);
    Field2d micro_field(GRID_N_MICRO);
    
    gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    int seed = time(NULL);
    gsl_rng_set(rng, seed);
    
    //std::cout << seed << std::endl;
    double macro_var = macro_field.expected_var(&macro_pl);
    double micro_var = micro_field.expected_var(&micro_pl);

    double total_var = macro_var + micro_var;
    //double total_var = macro_var;
    //double total_var = micro_var;
    
    CDF *gaussian = new CDF_GAUSSIAN(sqrt(total_var));

    std::cerr << "Macro field variance: " << macro_var << std::endl;
    std::cerr << "Micro field variance: " << micro_var << std::endl;
    std::cerr << "Total field variance: " << total_var << std::endl;

    macro_field.fill_power_l(&macro_pl, rng);
    macro_field.fft_l_to_x();

    std::cerr << "Macro field actual variance: " << macro_field.var() << std::endl;

    for(int i = 0; i < GRID_N_MACRO; i++)
    {
        for(int j = 0; j < GRID_N_MACRO; j++)
        {
            micro_field.fill_power_l(&micro_pl, rng);
            micro_field.fft_l_to_x();
            //std::cerr << "Micro field initial variance: " << micro_field.var() << std::endl;
            micro_field.add(macro_field.real(i,j),macro_field.imag(i,j));
            micro_field.remap_cdf(gaussian, cdf);

            //std::cerr << "Micro field remapped variance: " << micro_field.var() << std::endl;

            double mean = micro_field.mean();

            std::cout << macro_field.real(i,j) << " " << mean << std::endl;

            macro_field.set(i,j,mean, 0.0);
            //std::cout << mean << std::endl;
            //std::cerr << mean << std::endl;
        }
    }

    //macro_field.remap_cdf(gaussian, cdf);

    std::cerr << "Macro field mean: " << macro_field.mean() << std::endl;
    std::cerr << "Macro field var: " << macro_field.var() << std::endl;

    /*std::cerr << "Constructed field variance: " << field.var() << std::endl;
    std::cerr << "Constructed field mean: " << field.mean() << std::endl;
    
    double var = field.expected_var(&p2_l, rng);*/
    
    //field.print_ps(&print);

    

    /*field.output_remap_and_reduce(gaussian, cdf, 4);

    field.remap_field_cdf(gaussian, cdf);
    
    std::cerr << "Remapped field variance: " << field.var() << std::endl;
    std::cerr << "Remapped field mean: " << field.mean() << std::endl;
    
    grid reduced_field(field, 4);

    std::cerr << "Reduced field variance: " << reduced_field.var() << std::endl;
    std::cerr << "Reduced field mean: " << reduced_field.mean() << std::endl;*/

    /*for(int i = 0; i < GRID_N; i++)
    {
        for(int j = 0; j < GRID_N; j++)
        {
            double val = field.get_real(i,j);
            //double real2 = shear2.get_real(i,j);

            //double val = sqrt(real1*real1+real2*real2);

            std::cout << val << std::endl;
        }
    }*/

    macro_field.fft_x_to_l();
    
    Field2d shear1(&macro_field);
    Field2d shear2(&macro_field);

    double e_crit = 118333.0;

    std::cerr << "log(e_crit): " << log(e_crit) << std::endl;

    shear1.apply_filter_shear1(e_crit);
    shear2.apply_filter_shear2(e_crit);

    shear1.fft_l_to_x();
    shear2.fft_l_to_x();

    std::cerr << "shear1 field variance: " << shear1.var() << std::endl;
    std::cerr << "shear1 field mean: " << shear1.mean() << std::endl;
    std::cerr << "shear2 field variance: " << shear2.var() << std::endl;
    std::cerr << "shear2 field mean: " << shear2.mean() << std::endl;

    /*for(int i = 0; i < GRID_N; i++)
    {
        for(int j = 0; j < GRID_N; j++)
        {
            double real1 = shear1.get_real(i,j);
            double real2 = shear2.get_real(i,j);

            double val = sqrt(real1*real1+real2*real2);

            std::cout << val << std::endl;
        }
    }*/

    double noise_sigma = 0.3 / sqrt(2.0) / sqrt(254.776) / (GRID_L_MACRO / GRID_N_MACRO) / sqrt(1000.0);

    std::cerr << "noise sigma: " << noise_sigma << std::endl;

    shear1.add_noise(noise_sigma, rng);
    shear2.add_noise(noise_sigma, rng);

    std::cerr << "after noise" << std::endl;

    std::cerr << "shear1 field variance: " << shear1.var() << std::endl;
    std::cerr << "shear1 field mean: " << shear1.mean() << std::endl;
    std::cerr << "shear2 field variance: " << shear2.var() << std::endl;
    std::cerr << "shear2 field mean: " << shear2.mean() << std::endl;

    shear1.fft_x_to_l();
    shear2.fft_x_to_l();

    shear1.reverse_filter_shear1(e_crit);
    shear2.reverse_filter_shear2(e_crit);

    shear1.fft_l_to_x();
    shear2.fft_l_to_x();

    std::cerr << "after reverse filtering" << std::endl;

    std::cerr << "shear1 field variance: " << shear1.var() << std::endl;
    std::cerr << "shear1 field mean: " << shear1.mean() << std::endl;
    std::cerr << "shear2 field variance: " << shear2.var() << std::endl;
    std::cerr << "shear2 field mean: " << shear2.mean() << std::endl;

    /*for(int i = 0; i < GRID_N_MACRO; i++)
    {
        for(int j = 0; j < GRID_N_MACRO; j++)
        {
            double real1 = shear1.real(i,j);
            double real2 = shear2.real(i,j);

            double val = sqrt(real1*real1+real2*real2);

            std::cout << val << std::endl;
        }
    }*/

    return 0;
}
