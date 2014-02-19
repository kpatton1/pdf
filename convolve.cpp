
#include <cmath>

#include <iostream>
#include <fstream>

#include "gsl/gsl_interp.h"
#include "gsl/gsl_integration.h"

const double sigma = 7.878 * 1.686;

struct param_struct
{
    int N;
    gsl_interp* interp;
    gsl_interp_accel* accel;
    
    double* E_p_array;
    double* prob_array;
    
    double E_p;
    
    double E_p_min;
    double E_p_max;
};

double gaussian(double x)
{
    double result = 1.0 / (sigma * sqrt(2.0 * M_PI)) * exp( -0.5 * x * x / (sigma * sigma));
    return result;
};

double integral(double x, void* params);

int main()
{
    std::ifstream ifs;
    ifs.open("halo_mod.dat", std::ifstream::in);
    
    double E_p_array[1000];
    double prob_array[1000];
    
    int i = 0;
    
    while(ifs.good() && i < 1000)
    {
        double E_p;
        double prob;
        ifs >> E_p >> prob;
        
        if(!ifs.good())
        {
            break;
        }
        
        if(i > 0 && E_p <= E_p_array[i-1])
        {
            continue;
        }
        
        E_p_array[i] = E_p;
        prob_array[i] = prob;
        
        i++;
        
    }
    
    std::cout << i << std::endl;
    
    gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, i);
    
    double E_p_min = E_p_array[0];
    double E_p_max = E_p_array[i-1];
    
    gsl_interp_init(interp, E_p_array, prob_array, i);
    
    gsl_interp_accel* accel = gsl_interp_accel_alloc();
    
    double integral_value = gsl_interp_eval_integ(interp, E_p_array, prob_array, E_p_min, E_p_max, accel);
    
    std::cout << integral_value << std::endl;

    for(int j = 0; j < i; j++)
    {
        prob_array[i] = prob_array[i] / integral_value;
    }
    
    gsl_integration_workspace* w;
    
    w = gsl_integration_workspace_alloc(10000);
    
    param_struct params;
    params.N = i;
    params.interp = interp;
    params.accel = accel;
    params.E_p_array = E_p_array;
    params.prob_array = prob_array;
    
    params.E_p_min = E_p_min;
    params.E_p_max = E_p_max;
    
    gsl_function f;
    f.function = &integral;
        
    f.params = &params;
    
    /*for(int j = 0; j < i; j++)
    {
        double result;
        double error;
        
        double E_p = E_p_array[j];
        
        params.E_p = E_p;
        
        gsl_integration_qag(&f, -3.0 * sigma, 3.0 * sigma, 1e-8, 1e-8, 10000, GSL_INTEG_GAUSS51, w, &result, &error);
        
        double prob = result * 0.02893 + (1.0 - 0.02893) * gaussian(E_p - 7.878);
        
        std::cout << E_p << " " << prob << std::endl;
        
    }*/

    gsl_integration_workspace_free(w);

    gsl_interp_accel_free(accel);
    gsl_interp_free(interp);
    
    
    return 0;
}

double integral(double x, void* params)
{
    param_struct* p = (param_struct*) params;
    double E_p = p->E_p;
    
    
    
    if(E_p + x < p->E_p_min || E_p + x > p->E_p_max)
    {
        return 0.0;
    }
    
    double pdf_prob = gsl_interp_eval(p->interp, p->E_p_array, p->prob_array, E_p + x, p->accel);
    double gauss_prob = gaussian(x);
    
    return pdf_prob * gauss_prob;
}

