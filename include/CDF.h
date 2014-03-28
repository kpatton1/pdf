#ifndef CDF_H
#define CDF_H

#include <iostream>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_interp.h"
#include <cmath>

class CDF
{
public:
    virtual double cdf(double lnx) { return 0.0; }
    
    virtual double cdf_inv(double cdf) { return 0.0; }
    
    virtual ~CDF() { }

};

class CDF_GAUSSIAN : public CDF
{
    double sigma;
public:
    CDF_GAUSSIAN(double sigma)
    {
        this->sigma = sigma;
    };

    double cdf(double x)
    {
        gsl_cdf_gaussian_P(x, this->sigma);
    };
    
    double cdf_inv(double cdf)
    {
    	return 0.0;
    };
};

class CDF_PDF_INTERP : public CDF
{
    int N;
    
    double* cdf_array;
    double* pdf_array;
    double* lnx_array;
    
    double lnx_min;
    double lnx_max;
    
    gsl_interp* pdf_interp;
    gsl_interp_accel* pdf_interp_accel;
  
    gsl_interp* cdf_interp;
    gsl_interp_accel* cdf_interp_accel;
  
    gsl_interp* cdf_inv_interp;
    gsl_interp_accel* cdf_inv_interp_accel;

public:
    CDF_PDF_INTERP(double* lnx_array, double* pdf_array, int N)
    {
        this->N = N;
        
        if(N > 0)
        {
            this->lnx_array = new double[N];
            this->pdf_array = new double[N];
            this->cdf_array = new double[N];
            
            lnx_min = lnx_array[0];
            lnx_max = lnx_array[N-1];
            
            for(int i = 0; i < N; i++)
            {
                this->pdf_array[i] = pdf_array[i];
                this->lnx_array[i] = lnx_array[i];
            }
                
            this->pdf_interp = gsl_interp_alloc(gsl_interp_cspline, N);
            gsl_interp_init(this->pdf_interp, this->lnx_array, this->pdf_array, N);
            this->pdf_interp_accel = gsl_interp_accel_alloc();
    
            for(int i = 0; i < N; i++)
            {
                this->cdf_array[i] = gsl_interp_eval_integ(this->pdf_interp, this->lnx_array, this->pdf_array, this->lnx_min, this->lnx_array[i], this->pdf_interp_accel);
                //std::cout << this->cdf_array[i] << std::endl;
            }
            
            for(int i = 0; i < N; i++)
            {
                this->cdf_array[i] = this->cdf_array[i] / this->cdf_array[N-1];
            }
           
            //std::cout << "test1" << std::endl;
          
            this->cdf_interp = gsl_interp_alloc(gsl_interp_cspline, N);
            gsl_interp_init(this->cdf_interp, this->lnx_array, this->cdf_array, N);
            this->cdf_interp_accel = gsl_interp_accel_alloc();
            
            //            std::cout << "test2" << std::endl;
            
            this->cdf_inv_interp = gsl_interp_alloc(gsl_interp_cspline, N);
            gsl_interp_init(this->cdf_inv_interp, this->cdf_array, this->lnx_array, N);
            this->cdf_inv_interp_accel = gsl_interp_accel_alloc();
            
            //            std::cout << "test3" << std::endl;
        }
        else
        {
            lnx_min = 0.0;
            lnx_max = 0.0;
            
            this->pdf_array = 0;
            this->lnx_array = 0;
            this->cdf_array = 0;
            
            this->pdf_interp = 0;
            this->pdf_interp_accel = 0;
            this->cdf_inv_interp = 0;
            this->cdf_inv_interp_accel = 0;
        }
        
        
        
    
    }

    double cdf(double x)
    {
        if(x <= 0.0)
            return 0.0;
        
        double lnx = log(x);
        
        if(lnx <= lnx_min)
            return 0.0;
            
        if(lnx >= lnx_max)
            return 1.0;
        
        double cdf = gsl_interp_eval(this->cdf_interp, this->lnx_array, this->cdf_array, lnx, this->cdf_interp_accel);
        
        return cdf;
    };
    
    double cdf_inv(double cdf)
    {
        if(cdf <= 0.0)
            return lnx_min;
        
        if(cdf >= 1.0)
            return lnx_max;
        
        double lnx = gsl_interp_eval(this->cdf_inv_interp, this->cdf_array, this->lnx_array, cdf, this->cdf_inv_interp_accel);
        
        return exp(lnx)-1.0;
    };
    
    ~CDF_PDF_INTERP()
    {
    
    };

};

#endif

