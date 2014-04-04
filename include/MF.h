#ifndef MF_H
#define MF_H

#include <cmath>
#include "PS.h"

class MF
{
public:
    virtual double dn_dlnm(double lnm, double z) { return 0.0; }
    virtual double mdn_dlnm(double lnm, double z) { return 0.0; }
    virtual ~MF(){}
};

class MF_PS : public MF
{
private:
    PS* ps;
    static const double dc = 1.686;
public:
    MF_PS(PS* ps)
    {
        this->ps = ps;
        
    }
    
    double dn_dlnm(double lnm, double z) // M * dN/dM as a function of (M/rho)
    {
        //double r = pow(m_p * 3.0 / (4.0 * M_PI), 1.0/3.0);
        double r = exp(1.0/3.0 * (lnm + log(3.0 / (4.0 * M_PI))));
        double m = exp(lnm);
        double var = this->ps->sigma2_tophat(r, z);
        double r_dvar_dr = this->ps->sigma2_r_dtophat(r, z);
        double f = sqrt(2.0 / (M_PI*var)) * dc * exp(-0.5 * dc * dc / var);
        
        double result = -0.5 * f * r_dvar_dr / (3.0 * m * var);
        
        return result;
    }

    double mdn_dlnm(double lnm, double z) // M * dN/dM as a function of (M/rho)
    {
        //double r = pow(m_p * 3.0 / (4.0 * M_PI), 1.0/3.0);
        double r = exp(1.0/3.0 * (lnm + log(3.0 / (4.0 * M_PI))));
        double var = this->ps->sigma2_tophat(r, z);
        double r_dvar_dr = this->ps->sigma2_r_dtophat(r, z);
        double f = sqrt(2.0 / (M_PI*var)) * dc * exp(-0.5 * dc * dc / var);
        
        double result = -0.5 * f * r_dvar_dr / (3.0 * var);
        
        return result;
    }
};

#endif

