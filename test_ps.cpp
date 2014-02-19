
#include <iostream>

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
    
    double L = 10.0;
    int N = 1000;
    
    double kmin = 2.0 * M_PI / L * 0.1;
    double kmax = 2.0 * M_PI / L * N * 10.0;
    
    int steps = 10000;
    
    for(double k = kmin; k < kmax; k += (kmax-kmin)/steps)
    {
        double pk = ps->p(k, 0.0);
        
        std::cout << k << " " << pk << std::endl;
    }
    
    return 0;
}
