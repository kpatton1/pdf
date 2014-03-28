
#include <iostream>
#include <cmath>

#include "gsl/gsl_interp.h"

#include "NFW.h"
#include "PS.h"
#include "TF.h"
#include "GF.h"
#include "MF.h"
#include "Halo.h"
#include "PDF.h"

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
    
    //double var = ps->var(0.01/h, 0.0);
    //std::cout << "gaussian field variance: " << var << std::endl;
    
    MF* mf = new MF_PS(ps);
    
    std::cout.precision(8);
    
    /*for(double m = 1.0001e10; m < 1e16; m *= 1.1)
    {
        double m_p = m / p;
        double r = pow(m_p * 3.0 / (4.0 * M_PI), 1.0/3.0);
        double var = ps->var(r, 0.0);
        
        
        
        double dn_dlogm = mf->m_dn_dm(m_p, 0.0);
        //std::cout << "m_p: " << m_p << "  dndm: " << dndm  << "  m_p*dndm: " << m_p*dndm << std::endl;
        std::cout << m << "\t" << r << "\t" << var << "\t" << dn_dlogm << std::endl;
    }*/
    
    /*for(double k = 0.01; k < 10.0; k += 0.01)
    {
        double pk = ps->p(k, 0.0);
        std::cout << k << " " << pk << std::endl;
    }*/
    
    double var = ps->sigma2_tophat(7.878, 0.0);
    
    double m = p * 4.0 * M_PI * 7.878 * 7.878 * 7.878 / 3.0;
    
    std::cerr << m << std::endl;
    
    std::cerr << var << std::endl;
    
    Halo* halo = new Halo_NFW(1e-10, 10.0, 1000000);
    PDF_HALO* pdf = new PDF_HALO(mf, halo);
    
    double lnm_min = log(1.0e0 / p);
    double lnm_max = log(1.0e18 / p);
    
    double covering = pdf->covering(lnm_min, lnm_max, 0.0);
    
    std::cerr << covering << std::endl;
    std::cerr << covering * 7.878 << std::endl;
    
    for(double lne = log(1.0e-4); lne < log(1.0e6); lne += log(1.2))
    {
        double pdfval = pdf->pdf(lne, lnm_min, lnm_max, 0.0);
        std::cout << lne << " " << pdfval / covering << std::endl;
        std::cerr << lne << " " << pdfval / covering << std::endl;
        //std::cout << e_p << " " << pdfval / covering << " " << log(e_p) << " " << e_p * pdfval / covering << std::endl;
    }
    

    
    
    
    
}


/*int main()
{
    const int BINS = 100000;
    double f[BINS];
    double x[BINS];
    
    double xnew[BINS];
    double dxdf[BINS];
    
    double min = log(1e-6);
    double max = log(1e6);
    double diff = (max - min) / BINS;
    for(int i = 0; i < BINS; i++)
    {
        double xi = exp(max - diff*i);
        double fi = NFW_projected_mass(xi);

        f[i] = fi;        
        x[i] = xi;
    }
    
    gsl_interp* interp = gsl_interp_alloc(gsl_interp_cspline, BINS);
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_interp_init(interp, f, x, BINS);
    
    for(int i = 0; i < BINS; i++)
    {
        xnew[i] = gsl_interp_eval(interp, f, x, f[i], acc);
        dxdf[i] = gsl_interp_eval_deriv(interp, f, x, f[i], acc);
        
        std::cout << f[i] << " " << x[i] << " " << xnew[i] << " " << dxdf[i] << std::endl;
    }
    
    gsl_interp_accel_free(acc);
    gsl_interp_free(interp);
    
    return 0.;
}*/



