#ifndef HALO_H
#define HALO_H


#include "gsl/gsl_interp.h"
#include <cmath>
#include <iostream>

class Halo
{
public:
    virtual double PDF(double lne, double lnm, double c) { return 0.0; };
/*    double (*f)(double);
    
    int N;
    double* lnf_arr;
    double* x_arr;
    
    double xmin;
    double xmax;
    
    double fmin;
    double fmax;
    
    gsl_interp* interp;
    gsl_interp_accel* acc;
    

public:
    Halo()
    {
        this->f = 0;
        this->N = 0;
        this->xmin = 1.0;
        this->xmax = 1.0;
        this->fmin = 1.0;
        this->fmax = 1.0;
    
        this->lnf_arr = 0;
        this->x_arr = 0;
        this->interp = 0;
        this->acc = 0;
    }
    
    Halo( double (*f)(double), double xmin, double xmax, int N )
    {
        this->f = f;
        this->xmin = min;
        this->xmax = max;
        this->N = N;
        
        double logmin = log(min);
        double logmax = log(max);
        
        double logstep = (logmax - logmin) / N;
        
        this->lnf_arr = new double[N];
        this->x_arr = new double[N];
        
        this->interp = gsl_interp_alloc(gsl_interp_cspline, N);
        this->acc = gsl_interp_accel_alloc();
        
        if(!this->lnf_arr || !this->x_arr || !this->interp || !this->acc)
        {
            std::cout << "Error allocating memory for Halo interpolation!" << std::endl;
            this->N = 0;
            this->xmin = 1.0;
            this->xmax = 1.0;
            this->fmin = 1.0;
            this->fmax = 1.0;
            return;
        }
        
        for(int i = 0; i < this->N; i++)
        {
            double xi = exp(logmax - logstep * i);
            x_arr[i] = xi; 
            double lnfi = log(f(xi));
            lnf_arr[i] = lnfi;
            
            //std::cout << xi << " " << fi << std::endl;
        }
        
        this->fmin = lnf_arr[0];
        this->fmax = lnf_arr[N-1];
        
        //std::cout << "halo init fmin: " << fmin << "  fmax: " << fmax << std::endl;
        
        gsl_interp_init(this->interp, this->lnf_arr, this->x_arr, N);
    }
    
    virtual double PDF(double lne, double lnm, double c)
    {
        double em_p = 2.0 * pow(m_p * 3.0 / (800.0 * M_PI), 1.0/3.0) * 200.0 / 3.0 * c * c * c / (log(1.0 + c) - c / (1.0 + c));

        double f = e_p/em_p;
        
        if(f < fmin)
            return 0.0;
           
        if(f > fmax)
        {
            return 0.0;
        }
        
    //    std::cout << e_p << " " << m_p << " " << em_p << " " << f << std::endl;
        
        double x = gsl_interp_eval(this->interp, this->lnf_arr, this->x_arr, f, this->acc);
        double dxdf = gsl_interp_eval_deriv(this->interp, this->lnf_arr, this->x_arr, f, this->acc);
        
    //    double xcut = gsl_interp_eval(this->interp_xf, this->f_arr, this->x_arr, f, this->acc_xf);
        double xcut = c;
        
        double pdf = (-2.0) * x * dxdf / (em_p * c * c);
        
        return pdf;
    }

    ~Halo()
    {
        if(this->acc)
        {
            gsl_interp_accel_free(this->acc);
        }
        if(this->interp)
        {
            gsl_interp_free(this->interp);
        }
        if(this->lnf_arr)
        {
            free(this->lnf_arr);
        }
        if(this->x_arr)
        {
            free(this->x_arr);
        }
    }
*/
};

class Halo_NFW : public Halo
{
    
    double f(double x)
    {
        if(x == 1.0)
        {
            return 1.0/3.0;
        }
        if(x < 1.0)
        {   
           double ret = 1.0 / (x*x - 1.0) * (1.0 - 2.0 / sqrt(1.0 - x*x) * atanh(sqrt((1.0-x)/(x+1.0))));
           
           return ret;
        }
        else
        {
           double ret = 1.0 / (x*x - 1.0) * (1.0 - 2.0 / sqrt(x*x-1.0) * atan(sqrt((x-1.0)/(x+1.0))));
           
           return ret;
        }
    }
    
    int N;
    double* lnf_arr;
    double* x_arr;
    
    double xmin;
    double xmax;
    
    double lnfmin;
    double lnfmax;
    
    gsl_interp* interp;
    gsl_interp_accel* acc;
    
public:
    Halo_NFW(double xmin, double xmax, int N)
    {
        this->xmin = xmin;
        this->xmax = xmax;
        this->N = N;
        
        double logmin = log(xmin);
        double logmax = log(xmax);
        
        double logstep = (logmax - logmin) / N;
        
        this->lnf_arr = new double[N];
        this->x_arr = new double[N];
        
        this->interp = gsl_interp_alloc(gsl_interp_cspline, N);
        this->acc = gsl_interp_accel_alloc();
        
        if(!this->lnf_arr || !this->x_arr || !this->interp || !this->acc)
        {
            std::cout << "Error allocating memory for Halo interpolation!" << std::endl;
            this->N = 0;
            this->xmin = 1.0;
            this->xmax = 1.0;
            this->lnfmin = 0.0;
            this->lnfmax = 0.0;
            return;
        }
        
        for(int i = 0; i < this->N; i++)
        {
            double xi = exp(logmax - logstep * i);
            x_arr[i] = xi; 
            double lnfi = log(f(xi));
            lnf_arr[i] = lnfi;
            
            //std::cout << xi << " " << fi << std::endl;
        }
        
        this->lnfmin = lnf_arr[0];
        this->lnfmax = lnf_arr[N-1];
        
        //std::cout << "halo init fmin: " << fmin << "  fmax: " << fmax << std::endl;
        
        gsl_interp_init(this->interp, this->lnf_arr, this->x_arr, N);     
    }
    
    ~Halo_NFW()
    {
        if(this->acc)
        {
            gsl_interp_accel_free(this->acc);
        }
        if(this->interp)
        {
            gsl_interp_free(this->interp);
        }
        if(this->lnf_arr)
        {
            free(this->lnf_arr);
        }
        if(this->x_arr)
        {
            free(this->x_arr);
        }
    }
    
    double PDF(double lne, double lnm, double c)
    {
        // double em_p = 2.0 * pow(m_p * 3.0 / (800.0 * M_PI), 1.0/3.0) * 200.0 / 3.0 * c * c * c / (log(1.0 + c) - c / (1.0 + c));
        double lnem = 1.0/3.0 * (lnm + log(3.0 / (800.0 * M_PI))) + log (2.0 * 200.0 / 3.0) + log(c*c / (log(1.0 + c) - c / (1.0 + c)));

        double lnf = lne - lnem;
    
        if(lnf < lnfmin)
            return 0.0;
        
        if(lnf > lnfmax)
        {
            return 0.0;
        }
        
    //    std::cout << e_p << " " << m_p << " " << em_p << " " << f << std::endl;
        
        double x = gsl_interp_eval(this->interp, this->lnf_arr, this->x_arr, lnf, this->acc);
        double dxdf = gsl_interp_eval_deriv(this->interp, this->lnf_arr, this->x_arr, lnf, this->acc);
        
    //    double xcut = gsl_interp_eval(this->interp_xf, this->f_arr, this->x_arr, f, this->acc_xf);
        double xcut = c;
        
        double pdf = (-2.0) * x * dxdf / (c * c);
        
        return pdf;
    }
    
};

#endif

