#ifndef PS_H
#define PS_H

#include <cmath>
#include "gsl/gsl_integration.h"

#include "TF.h"
#include "GF.h"


double filter_tophat_k(double k, double r);
double dfilter_tophat_k_dr(double k, double r);

double filter_tophat_u(double u);
double dfilter_tophat_u_du(double u);

double gsl_integral_var(double k, void* params);
double gsl_integral_dvar_dr(double k, void* params);

double gsl_integral_var_rescaled(double x, void* params);
double gsl_integral_dvar_dr_rescaled(double x, void* params);

double gsl_norm_integral_ps_lss(double k, void* params);

double gsl_norm_integral_ps_lss_rescaled(double u, void* params);

struct var_struct
{
    void* ps;
    double r;
    double z;
};

class PS
{
protected:
    //gsl_integration_workspace* w;
    gsl_integration_cquad_workspace* w;
    
public:
    PS()
    {
        //this->w = gsl_integration_workspace_alloc(1000000);
        this->w = gsl_integration_cquad_workspace_alloc(1000);
    };

    virtual double p(double k, double z) { return 0.0; }; // delta^2 (k) - fractional variance per log k
    
    double var(double r, double z)
    {
        double result;
        double error;
        size_t nevals;
        
        gsl_function f;
        //f.function = &gsl_integral_var;
        f.function = &gsl_integral_var_rescaled;
        
        var_struct params;
        params.ps = this;
        params.r = r;
        params.z = z;
        
        f.params = &params;
        
        //gsl_integration_qagiu(&f, 0.0, 1e-8, 1e-6, 1000000, w, &result, &error);
        gsl_integration_cquad(&f, 0.0, 1.0, 1e-9, 1e-7, this->w, &result, &error, &nevals);
        
        return result;
    }
    
    double dvar_dr(double r, double z)
    {
        double result;
        double error;
        size_t nevals;
        
        gsl_function f;
        //f.function = &gsl_integral_dvar_dr;
        f.function = &gsl_integral_dvar_dr_rescaled;
        
        var_struct params;
        params.ps = this;
        params.r = r;
        params.z = z;
        
        f.params = &params;
        
        //gsl_integration_qagiu(&f, 0.0, 1e-8, 1e-6, 1000000, w, &result, &error);
        gsl_integration_cquad(&f, 0.0, 1.0, 1e-9, 1e-7, this->w, &result, &error, &nevals);
        
        return result/r;
    }
    
    ~PS()
    {
        if(this->w)
        {
            //gsl_integration_workspace_free(this->w);
            gsl_integration_cquad_workspace_free(this->w);
            this->w = 0;
        }
    };
};

class PS_LSS : public PS
{
private:
    TF* tf;
    GF* gf;
    double s8;
    double ns;
    double h;
    
    double ps_norm;

public:
    friend double gsl_norm_integral_ps_lss(double k, void* params);
    friend double gsl_norm_integral_ps_lss_rescaled(double u, void* params);
    
    PS_LSS(TF* tf, GF* gf, double h, double s8, double ns)
    {
        this->tf = tf;
        this->gf = gf;
        this->s8 = s8;
        this->ns = ns;
        this->h = h;
        
        double result;
        double error;
        size_t nevals;
        
        gsl_function f;
        f.function = &gsl_norm_integral_ps_lss_rescaled;
        f.params = this;
        
        //gsl_integration_qagiu(&f, 0.0, 1e-8, 1e-8, 10000, w, &result, &error);
        gsl_integration_cquad(&f, 0.0, 1.0, 1e-9, 1e-7, this->w, &result, &error, &nevals);
        
        this->ps_norm = s8 * s8 / result;
        
//        std::cout << "ret: " << ret << "  result: " << result << "  error: " << error << std::endl;
//        std::cout << "ps_norm: " << this->ps_norm << std::endl;
    };
    
    double p(double k, double z)
    {
        double tk_k2 = this->tf->t_k2(k);
        
        if (tk_k2 * tk_k2 == 0.0)
        {
            return 0.0;
        }
        
        double result = tk_k2 * tk_k2 * this->ps_norm * pow(k, this->ns - 1.0);
        
//        std::cout << "eval p(k) k=" << k << "  ps_norm=" << ps_norm << "  tk=" << tk << "  result=" << result << std::endl;
        return result;
    };
    
    ~PS_LSS()
    {

    };
};

double filter_tophat_k(double k, double r)
{
    double y = k * r;
    
    double result = 3.0 * (sin(y) - y * cos(y)) / (y * y * y);
    
    return result;
}



double dfilter_tophat_k_dr(double k, double r)
{ 
    double y = k * r;
    double y2 = y * y;
    double y4 = y2 * y2;
    
    double result = 3.0 * k * (y2 * sin(y) + 3.0 * y * cos(y) - 3.0 * sin(y)) / y4;
    
    return result;
}

double filter_tophat_u(double u)
{
    double result = 3.0 * (sin(u) - u * cos(u)) / (u * u * u);
    
    return result;
}

double dfilter_tophat_u_du(double u)
{
    double u2 = u * u;
    double u4 = u2 * u2;
    
    double result = 3.0 * (u2 * sin(u) + 3.0 * u * cos(u) - 3.0 * sin(u)) / u4;
    
    return result;
}

double gsl_integral_var(double k, void* params)
{
    var_struct* p = (var_struct*) params;
    PS* ps = (PS*) p->ps;
    double r = p->r;
    double z = p->z;
    double pk = ps->p(k, z);
    double fk = filter_tophat_k(k, r);
    
    double result = fk * fk * pk / k;
    
//    std::cout << "eval int k=" << k << "  pk=" << pk << "  result=" << result << std::endl;
    
    return result;
}

double gsl_integral_dvar_dr(double k, void* params)
{
    var_struct* p = (var_struct*) params;
    PS* ps = (PS*) p->ps;
    double r = p->r;
    double z = p->z;
    double pk = ps->p(k, z);
    double fk = filter_tophat_k(k, r);
    double dfk_dr = dfilter_tophat_k_dr(k, r);
    
    double result = 2.0 * dfk_dr * fk * pk / k;
    
//    std::cout << "eval int k=" << k << "  pk=" << pk << "  result=" << result << std::endl;
    
    return result;
}

double gsl_integral_var_rescaled(double x, void* params)
{
    var_struct* p = (var_struct*) params;
    PS* ps = (PS*) p->ps;
    
    double u = x / (1.0 - x);
    
    double r = p->r;
    double z = p->z;
    double pk = ps->p(u/r, z);
    double fk = filter_tophat_u(u);
    
    double result = fk * fk * pk / u / (1.0- x) / (1.0-x);
    
//    std::cout << "eval int k=" << k << "  pk=" << pk << "  result=" << result << std::endl;
    
    return result;
}

double gsl_integral_dvar_dr_rescaled(double x, void* params)
{
    var_struct* p = (var_struct*) params;
    PS* ps = (PS*) p->ps;
    
    double u = x / (1.0 - x);
    
    double r = p->r;
    double z = p->z;
    double pk = ps->p(u/r, z);
    double fk = filter_tophat_u(u);
    double dfk_du = dfilter_tophat_u_du(u);
    
    double result = 2.0 * dfk_du * fk * pk / (1.0- x) / (1.0-x);
    
//    std::cout << "eval int k=" << k << "  pk=" << pk << "  result=" << result << std::endl;
    
    return result;
}



double gsl_norm_integral_ps_lss(double k, void* params)
{
    PS_LSS* ps = (PS_LSS*) params;
    double expk = pow(k, ps->ns-1.0) / k; // pow(k, ps->ns+3.0) / k
    double fk = filter_tophat_k(k, 8.0 / ps->h);
    double tk_k2 = ps->tf->t_k2(k);
    
    double result = expk * fk * fk * tk_k2 * tk_k2;
    
//    std::cout << "eval k=" << k << "  fk=" << fk << "  tk=" << tk << "  result=" << result << std::endl;
    
    return result;
}

double gsl_norm_integral_ps_lss_rescaled(double u, void* params)
{
    PS_LSS* ps = (PS_LSS*) params;
    
    double k = u / (1.0 - u);
    
    double expk = pow(k, ps->ns-1.0) / k; // pow(k, ps->ns+3.0) / k
    double fk = filter_tophat_k(k, 8.0 / ps->h);
    double tk_k2 = ps->tf->t_k2(k);
    
    double result = expk * fk * fk * tk_k2 * tk_k2 / (1.0- u) / (1.0-u);
    
//    std::cout << "eval k=" << k << "  fk=" << fk << "  tk=" << tk << "  result=" << result << std::endl;
    
    return result;
}



#endif

