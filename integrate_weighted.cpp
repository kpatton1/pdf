#include <cmath>

#include <iostream>
#include <fstream>

#include "gsl/gsl_interp.h"

int main()
{
    double x_array[100000];
    double y_array[100000];
    
    int i = 0;
    
    while(std::cin.good() && i < 100000)
    {
        double x;
        double y;
        std::cin >> x >> y;
        
        if(!std::cin.good())
        {
            break;
        }
        
        if(i > 0 && x <= x_array[i-1])
        {
            continue;
        }
        
        x_array[i] = x;
        y_array[i] = y * exp(x);
        
        i++;
        
    }
    
    //std::cout << i << std::endl;
    
    gsl_interp* interp = gsl_interp_alloc(gsl_interp_cspline, i);
    
    double x_min = x_array[0];
    double x_max = x_array[i-1];
    
    gsl_interp_init(interp, x_array, y_array, i);
    
    gsl_interp_accel* accel = gsl_interp_accel_alloc();
  
    //std::cout << integral_value << std::endl;
    
    for(int j = 0; j < i; j++)
    {
        double cdf = gsl_interp_eval_integ(interp, x_array, y_array, x_min, x_array[j], accel);
        std::cout << x_array[j] << " " << cdf << std::endl;
//        std::cout << cdf << " " << exp(E_p_array[j]) << std::endl;
    }

    gsl_interp_accel_free(accel);
    gsl_interp_free(interp);
    
    
    return 0;
}
