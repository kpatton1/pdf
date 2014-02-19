
#include "gsl/gsl_fft_complex.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_interp.h"

#include "CDF.h"

#include <iostream>

class grid
{
private:

    int N;

    gsl_complex* data;
    
public:

    grid(int N)
    {
        this->data = new gsl_complex[N*N];
        this->N = N;
    };

    ~grid()
    {
        if(this->data)
        {
            free(this->data);
            this->data = NULL;
        }
    };

    void fft_l_to_x()
    {
        gsl_fft_complex_wavetable* wavetable;
        gsl_fft_complex_workspace* workspace;
        
        wavetable = gsl_fft_complex_wavetable_alloc(this->N);
        workspace = gsl_fft_complex_workspace_alloc(this->N);
        
        for(int i = 0; i < this->N; i++)
        {
            gsl_fft_complex_forward((double*) &this->data[i*this->N], 1, this->N, wavetable, workspace);
        }
        
        for(int j = 0; j < this->N; j++)
        {
            gsl_fft_complex_forward((double*) &this->data[j], this->N, this->N, wavetable, workspace);
        }
        
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
        
    };
    
    
    void fft_x_to_l()
    {
        gsl_fft_complex_wavetable* wavetable;
        gsl_fft_complex_workspace* workspace;
        
        wavetable = gsl_fft_complex_wavetable_alloc(this->N);
        workspace = gsl_fft_complex_workspace_alloc(this->N);
        
        for(int i = 0; i < this->N; i++)
        {
            gsl_fft_complex_inverse((double*) &this->data[i*this->N], 1, this->N, wavetable, workspace);
        }
        
        for(int j = 0; j < this->N; j++)
        {
            gsl_fft_complex_inverse((double*) &this->data[j], this->N, this->N, wavetable, workspace);
        }
        
        gsl_fft_complex_wavetable_free(wavetable);
        gsl_fft_complex_workspace_free(workspace);
        
    };

    void set(int i, int j, double val1, double val2)
    {
        GSL_SET_COMPLEX(&this->data[i*this->N+j], val1, val2);
    };

/*    void fill( double (*f)(double) )
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real, imag;
                
                double x, y;
                
                //x = gsl_freq_map(i, this->N);
                //y = gsl_freq_map(j, this->N);
                x=i;
                y=j;
                
                double l = sqrt(x*x + y*y);
                
                double Pl = f(l);
                
                Pl = sqrt(Pl);
                
                real = Pl;
                imag = 0.0;
                
                this->set(i,j,real,imag);
            }
        }
    };*/
    
    void fill_random_l( double(*f)(double), gsl_rng* rng)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real, imag;
                
                double x, y;
                
                //x = gsl_freq_map(i, this->N);
                //y = gsl_freq_map(j, this->N);
                x=i;
                y=j;
                
                double l = sqrt(x*x + y*y);
                
                double Pl = f(l);
                
                Pl = sqrt(Pl);
                
                gsl_ran_dir_2d(rng, &x, &y);
                
                real = x*Pl;
                imag = y*Pl;
                
                this->set(i,j,real,imag);
            }
        }
    };
    
    void print_ps( void(*print)(double, double))
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real, imag;
                
                double x, y;
                
                //x = gsl_freq_map(i, this->N);
                //y = gsl_freq_map(j, this->N);
                x=i;
                y=j;
                
                double l = sqrt(x*x + y*y);
                
                real = GSL_REAL(this->data[i*this->N+j]);
                imag = GSL_IMAG(this->data[i*this->N+j]);
                
                double Pl = real * real + imag * imag;
                //double Pl = real * real;
                
                print(l, Pl);
            }
        }
    };
    
    void remap_field_cdf(CDF* gaussian, CDF* target)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double temp = GSL_REAL(this->data[i*this->N+j]);
                //double imag = GSL_IMAG(this->data[i*this->N+j]);
                
                //double temp = real + imag;
                
                double cdf = gaussian->cdf(temp);
                
                double x = target->cdf_inv(cdf);
                
                this->set(i,j,x,0.0);
            }
        }
    };
    
    double expected_var( double(*f)(double), gsl_rng* rng)
    {
        double xx = 0.0;
        
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real, imag;
                
                double x, y;
                
                //x = gsl_freq_map(i, this->N);
                //y = gsl_freq_map(j, this->N);
                x=i;
                y=j;
                
                double l = sqrt(x*x + y*y);
                
                double Pl = f(l);
                
                //Pl = sqrt(Pl);
                
                xx += Pl;
            }
        }
        
        return xx * 0.5;
    };
    
    
    double var()
    {
        double x = 0.0;
        double xx = 0.0;
        
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double temp = GSL_REAL(this->data[i*this->N+j]);
                //double imag = GSL_IMAG(this->data[i*this->N+j]);
                
                //std::cout << temp << std::endl;
                
                //double temp = sqrt(real*real + imag*imag);
                
                x+=temp;
                xx+=temp*temp;
            }
        }
        
        x = x / this->N / this->N;
        
        xx = xx / this->N / this->N;
        
        xx = xx - x*x;
        
        return xx;
    };

    double mean()
    {
        double x = 0.0;
        //double xx = 0.0;
        
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double temp = GSL_REAL(this->data[i*this->N+j]);
                
                //std::cout << temp << std::endl;
                
                x += temp;
                //xx += temp*temp;
            }
        }
        
        x = x / this->N / this->N;
        
        //xx = xx / this->N / this->N;
        
        //xx = xx - x*x;
        
        return x;
    };


};
