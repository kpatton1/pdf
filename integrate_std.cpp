#include <cmath>

#include <iostream>
#include <fstream>

#include "gsl/gsl_interp.h"

int main()
{
    double E_p_array[100000];
    double prob_array[100000];
    
    int i = 0;
    
    while(std::cin.good() && i < 100000)
    {
        double E_p;
        double prob;
        std::cin >> E_p >> prob;
        
        if(!std::cin.good())
        {
            break;
        }
        
        if(i > 0 && exp(E_p) <= E_p_array[i-1])
        {
            continue;
        }
        
        E_p_array[i] = exp(E_p);
        prob_array[i] = prob / (exp(E_p)) * prob / (exp(E_p));
        
        i++;
        
    }
    
    //std::cout << i << std::endl;
    
    gsl_interp* interp = gsl_interp_alloc(gsl_interp_cspline, i);
    
    double E_p_min = E_p_array[0];
    double E_p_max = E_p_array[i-1];
    
    gsl_interp_init(interp, E_p_array, prob_array, i);
    
    gsl_interp_accel* accel = gsl_interp_accel_alloc();
  
    //std::cout << integral_value << std::endl;
    
    for(int j = 0; j < i; j++)
    {
        double cdf = gsl_interp_eval_integ(interp, E_p_array, prob_array, E_p_min, E_p_array[j], accel);
        std::cout << E_p_array[j] << " " << cdf << std::endl;
//        std::cout << cdf << " " << exp(E_p_array[j]) << std::endl;
    }

    gsl_interp_accel_free(accel);
    gsl_interp_free(interp);
    
    
    return 0;
}
