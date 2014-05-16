
#include <iostream>
#include <cmath>

#include "PS.h"
#include "TF.h"
#include "GF.h"

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
    //PS* ps_nonlin = new PS_HALOFIT(ps, Om);
    
    double L = 10.0;
    int N = 1000;
    
    //double kmin = 2.0 * M_PI / L * 0.001;
    //double kmax = 2.0 * M_PI / L * N * 10.0;
    double logkmin = log(1/(1000.0));
    double logkmax = log(1/(0.0001));
    
    int steps = 1000;
    
    for(double logk = logkmin; logk < logkmax; logk += (logkmax-logkmin)/steps)
    {
        double k = exp(logk);
        double pk = ps->p(k, 0.0);
        //double pk_nonlin = ps_nonlin->p(k, 0.0);
        
        //std::cout << k << " " << pk << " " << pk_nonlin << std::endl;
        std::cout << k << " " << pk << std::endl;
    }
    
    return 0;
}
