#ifndef PS_H
#define PS_H

#include <cmath>
#include "gsl/gsl_integration.h"

#include "TF.h"
#include "GF.h"


struct filter
{
    double (*f)(double kr);
    double (*df)(double kr);
};

double filter_tophat_f(double kr)
{
    double result = 3.0 * (sin(kr) - kr * cos(kr)) / (kr * kr * kr);

    return result;
}

double filter_tophat_df(double kr)
{
    double kr2 = kr * kr;
    double kr4 = kr2 * kr2;

    double result = 3.0 * (kr2 * sin(kr) + 3.0 * kr * cos(kr) - 3.0 * sin(kr)) / kr4;

    return result;
}

filter filter_tophat = {filter_tophat_f, filter_tophat_df};

static double filter_gaussian_f(double kr)
{
    double kr2 = kr * kr;

    double result = exp(-kr2);

    return result;
}

static double filter_gaussian_df(double kr)
{
    double kr2 = kr * kr;

    double result = exp(-kr2) * (-2.0 * kr);

    return result;
}

filter filter_gaussian = {filter_gaussian_f, filter_gaussian_df};

class PS
{
protected:
    gsl_integration_cquad_workspace* w;



    struct sigma2_filtered_integral_struct
    {
        PS* ps;
        struct filter* filter;
        double r;
        double z;
    };

    static double sigma2_filtered_integral(double x, void* params)
    {
        sigma2_filtered_integral_struct* s = (sigma2_filtered_integral_struct*) params;

        double kr = x / (1.0 - x);

        double r = s->r;
        double z = s->z;

        double fk = s->filter->f(kr);
        double pk = s->ps->p(kr / r, z);

        double result = fk * fk * pk / kr / (1.0 - x) / (1.0 - x);

        return result;
    }

    static double dsigma2_filtered_integral(double x, void* params)
    {
        sigma2_filtered_integral_struct* s = (sigma2_filtered_integral_struct*) params;

        double kr = x / (1.0 - x);

        double r = s->r;
        double z = s->z;

        double fk = s->filter->f(kr);
        double dfk = s->filter->df(kr);
        double pk = s->ps->p(kr / r, z);

        double result = 2.0 * dfk * fk * pk / (1.0 - x) / (1.0 - x);

        return result;
    }





    double sigma2(double r, double z, filter* filter)
    {
        double result;
        double error;
        size_t nevals;

        gsl_function f;
        f.function = &sigma2_filtered_integral;

        sigma2_filtered_integral_struct params;
        params.ps = this;
        params.filter = filter;
        params.r = r;
        params.z = z;

        f.params = &params;

        gsl_integration_cquad(&f, 0.0, 1.0, 1e-9, 1e-7, this->w, &result, &error, &nevals);

        return result;
    }

    double dsigma2(double r, double z, filter* filter)
    {
        double result;
        double error;
        size_t nevals;

        gsl_function f;
        f.function = &dsigma2_filtered_integral;

        sigma2_filtered_integral_struct params;
        params.ps = this;
        params.filter = filter;
        params.r = r;
        params.z = z;

        f.params = &params;

        gsl_integration_cquad(&f, 0.0, 1.0, 1e-9, 1e-7, this->w, &result, &error, &nevals);

        return result / r;
    }

public:
    PS()
    {
        this->w = gsl_integration_cquad_workspace_alloc(1000);
    }

    virtual double p(double k, double z) // delta^2 (k) - fractional variance per log k
    {
        return 0.0;
    }

    double sigma2_tophat(double r, double z)
    {
        return sigma2(r,z,&filter_tophat);
    }

    double dsigma2_tophat(double r, double z)
    {
        return dsigma2(r,z,&filter_tophat);
    }

    double sigma2_gaussian(double r, double z)
    {
        return sigma2(r,z,&filter_gaussian);
    }

    double dsigma2_gaussian(double r, double z)
    {
        return dsigma2(r,z,&filter_gaussian);
    }

    virtual ~PS()
    {
        if (this->w)
        {
            gsl_integration_cquad_workspace_free(this->w);
            this->w = 0;
        }
    }
};

class PS_2D: public PS
{
private:
    PS* ps;
public:
    PS_2D(PS* ps)
    {
        this->ps = ps;
    }

    double p(double k, double z)
    {
        return 0.0;
    }
};

class PS_LSS: public PS
{
private:
    TF* tf;
    GF* gf;
    double s8;
    double ns;
    double h;

    double ps_norm;

    static double gsl_norm_integral_ps_lss(double k, void* params)
    {
        PS_LSS* ps = (PS_LSS*) params;
        double expk = pow(k, ps->ns - 1.0) / k; // pow(k, ps->ns+3.0) / k
        double fk = filter_tophat_f(k * 8.0 / ps->h);
        double tk_k2 = ps->tf->t_k2(k);

        double result = expk * fk * fk * tk_k2 * tk_k2;

        return result;
    }

    static double gsl_norm_integral_ps_lss_rescaled(double u, void* params)
    {
        PS_LSS* ps = (PS_LSS*) params;

        double k = u / (1.0 - u);

        double expk = pow(k, ps->ns - 1.0) / k; // pow(k, ps->ns+3.0) / k
        double fk = filter_tophat_f(k * 8.0 / ps->h);
        double tk_k2 = ps->tf->t_k2(k);

        double result = expk * fk * fk * tk_k2 * tk_k2 / (1.0 - u) / (1.0 - u);

        return result;
    }

public:

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
    }

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
    }

    ~PS_LSS()
    {

    }
};

class PS_HALOFIT: public PS
{
private:
    PS* ps;

public:
    PS_HALOFIT(PS* ps)
    {
        this->ps = ps;

    }

    double p(double k, double z)
    {
        return 0.0;
    }

    ~PS_HALOFIT()
    {

    }
};

double gsl_halofit_nonlinear(double lnk, void* params)
{

    return 0.0;
}

#endif

