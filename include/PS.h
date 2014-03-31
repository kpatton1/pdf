#ifndef PS_H
#define PS_H

#include <cmath>
#include "gsl/gsl_integration.h"
#include <gsl/gsl_roots.h>

#include "TF.h"
#include "GF.h"




class PS
{
protected:
    gsl_integration_cquad_workspace* w;

    static double filter_tophat_f(double kr)
    {
        double fk = 3.0 * (sin(kr) - kr * cos(kr)) / (kr * kr * kr);

        return fk*fk;
    }

    static double filter_tophat_df(double kr)
    {
        double fk = 3.0 * (sin(kr) - kr * cos(kr)) / (kr * kr * kr);

        double dfk = 3.0 * (kr*kr * sin(kr) + 3.0 * kr * cos(kr) - 3.0 * sin(kr)) / (kr*kr*kr*kr);

        return 2.0*fk*dfk*kr;
    }

    static double filter_gaussian_f(double kr)
    {
        double result = exp(-kr*kr);

        return result;
    }

    static double filter_gaussian_df(double kr)
    {
        double result = exp(-kr*kr) * (-2.0 * kr);

        return result*kr;
    }

    static double filter_gaussian_d2f(double kr)
    {
        double result = exp(-kr*kr) * (-2.0 - 4.0 * kr);

        return result*kr*kr;
    }

    struct sigma2_filtered_integral_struct
    {
        PS* ps;
        double (*filter)(double);
        double r;
        double z;
    };

    static double sigma2_filtered_integral(double x, void* params)
    {
        sigma2_filtered_integral_struct* s = (sigma2_filtered_integral_struct*) params;

        double kr = x / (1.0 - x);

        double r = s->r;
        double z = s->z;

        double fk = s->filter(kr);
        double pk = s->ps->p(kr / r, z);

        double result = fk * pk / kr / (1.0 - x) / (1.0 - x);

        return result;
    }

    double sigma2(double r, double z, double (*filter)(double))
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
        return sigma2(r,z,&filter_tophat_f);
    }

    double dsigma2_tophat(double r, double z)
    {
        return sigma2(r,z,&filter_tophat_df)/r;
    }

    double sigma2_gaussian(double r, double z)
    {
        return sigma2(r,z,&filter_gaussian_f);
    }

    double dsigma2_gaussian(double r, double z)
    {
        return sigma2(r,z,&filter_gaussian_df)/r;
    }

    double d2sigma2_gaussian(double r, double z)
    {
        return sigma2(r,z,&filter_gaussian_d2f)/r/r;
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

    double r_sigma = 8.0;

    double n_eff = 0.0;
    double C = 0.0;

    struct nonlinear_scale_struct
    {
        PS* ps;
        double z;
    };

    static double nonlinear_scale(double r, void* params)
    {
        nonlinear_scale_struct* s = (nonlinear_scale_struct*) params;

        return s->ps->sigma2_gaussian(r, s->z) - 1.0;
    }

    const int max_iterations = 1000;

public:
    PS_HALOFIT(PS* ps, double z)
    {
        this->ps = ps;

        gsl_function F;
        nonlinear_scale_struct s;

        F.function = nonlinear_scale;
        s.ps = this->ps;
        s.z = z;

        F.params = &s;

        gsl_root_fsolver *solver;
        solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        gsl_root_fsolver_set(solver, &F, 1.0e-3, 1.0e3);

        int iter = 0;
        while(iter < max_iterations)
        {
            int status = gsl_root_fsolver_iterate(solver);

            if (status == GSL_EBADFUNC)
            {
                std::cerr << "Error finding nonlinear scale for halofit: GSL_EBADFUNC" << std::endl;
                std::cerr << "at value r: " << gsl_root_fsolver_root(solver) << std::endl;
                break;
            }

            if (status == GSL_EBADFUNC)
            {
                std::cerr << "Error finding nonlinear scale for halofit: GSL_EZERODIV" << std::endl;
                std::cerr << "at value r: " << gsl_root_fsolver_root(solver) << std::endl;
                break;
            }

            double root = gsl_root_fsolver_root(solver);
            double val = nonlinear_scale(root, &s);

            status = gsl_root_test_residual(val, 1e-3);
            if(status != GSL_CONTINUE)
            {
                r_sigma = val;
                break;
            }

            iter++;
            if(iter > max_iterations)
            {
                std::cerr << "Error finding nonlinear scale for halofit: did not converge after " << iter << " iterations" << std::endl;
                std::cerr << "at value r: " << root << "  f(r): " << val << std::endl;
                break;
            }
        }

        gsl_root_fsolver_free(solver);

    }

    double p(double k, double z)
    {
        double p_lin = ps->p(k,z);

        double y = k * r_sigma;

        double fy = y/4.0 + y*y/8.0;

        double neff = ps->dsigma2_gaussian();

        double omega_m;
        double omega_w;
        double w;

        double C;

        double f1;
        double f2;
        double f3;

        double a_n;
        double b_n;
        double c_n;
        double gamma_n;
        double alpha_n;
        double beta_n;
        double mu_n;
        double nu_n;

        double p_q = p_lin * (pow(1.0+p_lin,beta_n)/(1.0+alpha_n*p_lin)) * exp(-fy);
        double p_h = ((a_n * pow(y,3*f1))/(1+b_n*pow(y,f2)+pow(c_n*f3*y,3.0-gamma_n)))/(1.0+mu_n/y + nu_n/y/y);
    }

    ~PS_HALOFIT()
    {

    }
};



#endif

