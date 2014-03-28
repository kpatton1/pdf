
#include <iostream>
#include <cmath>

#include "Halo.h"
#include "NFW.h"

const double G = 4.302e-9; // G in units of Mpc Msolar^-1 (km/s)^2

int main()
{
    double s8 = 1.0;
    double h = 0.7;
    double Om = 0.30;
    double ns = 1.0;
    
    double pc = 3.0 * 100.0 * 100.0 * h * h / (8.0 * M_PI * G);
    double p = Om * pc;
    
    double m_p_min = 1.0e0 / p;
    double m_p_max = 1.0e16 / p;    
    
    Halo* halo = new Halo_NFW(1e-10, 10.0, 1000000);
    
    double factor;
    
    std::cin >> factor;
    
    double m_p = factor / p;
    
    double r200 = exp(1.0/3.0 * (log(m_p) + log(3.0 / (800.0 * M_PI))));
    
//    double m_p = factor * factor * factor * (800.0 * M_PI) / 3.0;
    
    for(double e_p = 1.0e-4; e_p < 1e10; e_p = e_p * 1.01)
    {
        double prob = halo->PDF(log(e_p), log(m_p), 10.0);
        
        std::cout << log(e_p) << " " << prob << std::endl;
    }


}
