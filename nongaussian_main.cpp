
#include <math.h>
#include "grid.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_rng.h"

#include "CDF.h"

#include <ctime>


#include <iostream>
#include <fstream>

#define GRID_N 512
#define GRID_L 10.0

#define ALPHA -1.0

#define PDF_SIZE 99
#define PS_SIZE 10000

double k_array[PS_SIZE];
double ps_array[PS_SIZE];

gsl_interp* ps_interp = gsl_interp_alloc(gsl_interp_cspline, PS_SIZE);
gsl_interp_accel* ps_interp_accel = gsl_interp_accel_alloc();

double p2_l(double l)
{
    if( l <= 0.0)
    {
        return 0.0;
    }

    double k = 2.0 * M_PI * l / GRID_L;
    
    if(k < k_array[0])
    {
        k = k_array[0];
        l = GRID_L * k / (2.0 * M_PI);
    }
    if(k > k_array[PS_SIZE-1])
    {
        k = k_array[PS_SIZE-1];
        l = GRID_L * k / (2.0 * M_PI);
    }
    
    double pl = gsl_interp_eval(ps_interp, k_array, ps_array, k, ps_interp_accel);
    pl = pl / (2.0 * M_PI  * l * l);
   
    return pl;
}

void print(double l, double pl)
{
    double k = 2.0 * M_PI * l / GRID_L;
    double pk = (2.0 * M_PI * l * l) * pl;
    
    std::cout << k << " " << pk << std::endl;
}

int main(int argc, char* argv[])
{

    double lnx_array[PDF_SIZE];
    double pdf_array[PDF_SIZE];
    
    std::ifstream ps_file("ps.dat", std::ifstream::in);
    std::ifstream pdf_file("pdf.dat", std::ifstream::in);
    
    for(int i = 0; i < PS_SIZE; i++)
    {
        if(!ps_file.good())
        {
            std::cerr << "Error reading ps.dat!" << std::endl;
            break;
        }
        
        ps_file >> k_array[i] >> ps_array[i];  
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
    
    ps_interp = gsl_interp_alloc(gsl_interp_cspline, PS_SIZE);
    ps_interp_accel = gsl_interp_accel_alloc();
    gsl_interp_init(ps_interp, k_array, ps_array, PS_SIZE);
    
    std::cerr << "Set up ps interp!" << std::endl;
    
    CDF *cdf = new CDF_PDF_INTERP(lnx_array, pdf_array, PDF_SIZE);
    
    std::cout << "Set up CDF!" << std::endl;

    grid field(GRID_N);
    
    gsl_rng * rng = gsl_rng_alloc (gsl_rng_taus);
    int seed = time(NULL);
    gsl_rng_set(rng, seed);
    
    //std::cout << seed << std::endl;
    
    field.fill_random_l(&p2_l, rng);
    field.fft_l_to_x();
    
    std::cerr << "Expected gaussian variance: " << field.expected_var(&p2_l, rng) << std::endl;

    std::cerr << "Measured gaussian variance: " << field.var() << std::endl;
    std::cerr << "Measured gaussian mean: " << field.mean() << std::endl;
    
    double var = field.expected_var(&p2_l, rng);
    
    CDF *gaussian = new CDF_GAUSSIAN(sqrt(var));
    
    field.remap_field_cdf(gaussian, cdf);
    
    std::cerr << "Remapped variance: " << field.var() << std::endl;

    std::cerr << "Remapped mean: " << field.mean() << std::endl;
    
    field.fft_x_to_l();
    
    //field.print_ps(&print);

    return 0;
}
