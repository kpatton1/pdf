
#include <iostream>
#include <cmath>

#include "PS.h"
#include "TF.h"
#include "GF.h"
#include "MF.h"

const double G = 4.302e-9; // G in units of Mpc Msolar^-1 (km/s)^2

int main()
{
    double s8 = 1.0;
    double h = 0.7;
    double Om = 0.30;
    double ns = 1.0;
    
    double pc = 3.0 * 100.0 * 100.0 * h * h / (8.0 * M_PI * G);
    double p = Om * pc;
    
    //std::cerr << pc << std::endl;
    //std::cerr << p << std::endl;
    
    TF* tf = new TF_BBKS(Om, h);
    GF* gf = new GF_LCDM();
    
    PS* ps = new PS_LSS(tf, gf, h, s8, ns);
    
    MF* mf = new MF_PS(ps);
    
    double lnm_min = log(1.0e0 / p);
    double lnm_max = log(1.0e18 / p);
    
    for(double lnm = lnm_min; lnm < lnm_max; lnm += log(1.05))
    {
        double mdn_dlnm = mf->mdn_dlnm(lnm, 0.0);
        
        double r200 = exp(1.0/3.0 * (lnm + log(3.0 / (800.0 * M_PI))));
        
        std::cout << exp(lnm) * p << " " << (lnm+log(p))/log(10.0) << " " << log(10.0)*mdn_dlnm << " " << r200 << std::endl;
    }
    
    return 0;
}
