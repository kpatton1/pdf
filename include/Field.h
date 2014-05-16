#ifndef FIELD_H
#define FIELD_H

#include <cmath>
#include "gsl/gsl_fft_complex.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_interp.h"
#include "CDF.h"

class Field
{


};


class Field2d
{
private:
    int N;
    gsl_complex* data;

    gsl_fft_complex_wavetable* wavetable;
    gsl_fft_complex_workspace* workspace;

    double sinc_filter(double l)
    {
        if(l == 0.0)
            return 1.0;

        double x = M_PI * l / N;

        return sin(x) / x;
    }

public:
    Field2d(int N)
    {
        this->N = N;
        this->data = new gsl_complex[N*N];

        this->wavetable = gsl_fft_complex_wavetable_alloc(this->N);
        this->workspace = gsl_fft_complex_workspace_alloc(this->N);

        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                GSL_REAL(this->data[i*this->N+j]) = 0.0;
                GSL_IMAG(this->data[i*this->N+j]) = 0.0;
            }
        }
    }

    Field2d(Field2d *copy)
    {
        this->N = copy->N;
        this->data = new gsl_complex[N*N];

        this->wavetable = gsl_fft_complex_wavetable_alloc(this->N);
        this->workspace = gsl_fft_complex_workspace_alloc(this->N);

        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                GSL_REAL(this->data[i*this->N+j]) = GSL_REAL(copy->data[i*copy->N+j]);
                GSL_IMAG(this->data[i*this->N+j]) = GSL_IMAG(copy->data[i*copy->N+j]);
            }
        }
    }

    ~Field2d()
    {
        this->N = 0;

        if(this->data)
        {
            free(this->data);
            this->data = 0;
        }

        if(this->wavetable)
        {
            gsl_fft_complex_wavetable_free(wavetable);
            this->wavetable = 0;

        }

        if(this->workspace)
        {
            gsl_fft_complex_workspace_free(workspace);
            this->workspace = 0;
        }
    }

    double real(int i, int j)
    {
        return GSL_REAL(this->data[i*this->N+j]);
    }

    double imag(int i, int j)
    {
        return GSL_IMAG(this->data[i*this->N+j]);
    }

    void add(double real, double imag)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                GSL_REAL(this->data[i*this->N+j]) += real;
                GSL_IMAG(this->data[i*this->N+j]) += imag;
            }
        }
    }

    void set(int i, int j, double val1, double val2)
    {
        GSL_SET_COMPLEX(&this->data[i*this->N+j], val1, val2);
    }

    void fill_power_l( double (*P)(double l), gsl_rng* rng)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real, imag;

                double lx, ly;

                if(i > this->N / 2)
                {
                    lx = i - this->N;
                }
                else
                {
                    lx = i;
                }
                if(j > this->N / 2)
                {
                    ly = j - this->N;
                }
                else
                {
                    ly = j;
                }

                double l = sqrt(lx*lx + ly*ly);

                double Pl = P(l);

                Pl = sqrt(Pl) * sinc_filter(lx) * sinc_filter(ly);
                //Pl = sqrt(Pl);

                gsl_ran_dir_2d(rng, &real, &imag);

                real = real*Pl;
                imag = imag*Pl;

                this->set(i,j,real,imag);
            }
        }
    }

    void fft_x_to_l()
    {
        for(int i = 0; i < this->N; i++)
        {
            gsl_fft_complex_inverse((double*) &this->data[i*this->N], 1, this->N, wavetable, workspace);
        }

        for(int j = 0; j < this->N; j++)
        {
            gsl_fft_complex_inverse((double*) &this->data[j], this->N, this->N, wavetable, workspace);
        }
    }

    void fft_l_to_x()
    {
        for(int i = 0; i < this->N; i++)
        {
            gsl_fft_complex_forward((double*) &this->data[i*this->N], 1, this->N, wavetable, workspace);
        }

        for(int j = 0; j < this->N; j++)
        {
            gsl_fft_complex_forward((double*) &this->data[j], this->N, this->N, wavetable, workspace);
        }
    }

    void remap_cdf(CDF* gaussian, CDF* target)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double temp = GSL_REAL(this->data[i*this->N+j]);

                double cdf = gaussian->cdf(temp);

                double x = target->cdf_inv(cdf);

                this->set(i,j,x,0.0);
            }
        }
    }

    double expected_var( double(*f)(double))
    {
        double xx = 0.0;

        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real, imag;

                double lx, ly;

                if(i > this->N / 2)
                {
                    lx = i - this->N;
                }
                else
                {
                    lx = i;
                }
                if(j > this->N / 2)
                {
                    ly = j - this->N;
                }
                else
                {
                    ly = j;
                }

                double l = sqrt(lx*lx + ly*ly);

                double Pl = f(l) * sinc_filter(lx) * sinc_filter(ly) * sinc_filter(lx) * sinc_filter(ly);

                xx += Pl;
            }
        }

        return xx * 0.5;
    }


    double var()
    {
        double x = 0.0;
        double xx = 0.0;

        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double temp = GSL_REAL(this->data[i*this->N+j]);

                x+=temp;
                xx+=temp*temp;
            }
        }

        x = x / this->N / this->N;

        xx = xx / this->N / this->N;

        xx = xx - x*x;

        return xx;
    }

    double mean()
    {
        double x = 0.0;

        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double temp = GSL_REAL(this->data[i*this->N+j]);

                x += temp;
            }
        }

        x = x / this->N / this->N;

        return x;
    }

    void apply_filter_shear1(double e_crit)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real = GSL_REAL(this->data[i*this->N+j]);
                double imag = GSL_IMAG(this->data[i*this->N+j]);

                double lx, ly;

                if(i > this->N / 2)
                {
                    lx = i - this->N;
                }
                else
                {
                    lx = i;
                }
                if(j > this->N / 2)
                {
                    ly = j - this->N;
                }
                else
                {
                    ly = j;
                }

                if(i != 0 || j != 0)
                {
                    double filter = (lx*lx - ly*ly) / (lx*lx + ly*ly) / e_crit;
                    this->set(i,j,real*filter,imag*filter);
                }
                else
                {
                    double filter = 0.0;
                    this->set(i,j,real*filter,imag*filter);
                }
            }
        }
    }

    void reverse_filter_shear1(double e_crit)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real = GSL_REAL(this->data[i*this->N+j]);
                double imag = GSL_IMAG(this->data[i*this->N+j]);

                double lx, ly;

                if(i > this->N / 2)
                {
                    lx = i - this->N;
                }
                else
                {
                    lx = i;
                }
                if(j > this->N / 2)
                {
                    ly = j - this->N;
                }
                else
                {
                    ly = j;
                }

                if(i != 0 || j != 0)
                {
                    double filter = e_crit * (lx*lx - ly*ly) / (lx*lx + ly*ly);
                    //filter = filter * exp(-i*i - j*j);
                    this->set(i,j,real*filter,imag*filter);
                }
                else
                {
                    double filter = 0.0;
                    this->set(i,j,real*filter,imag*filter);
                }
            }
        }
    }

    void apply_filter_shear2(double e_crit)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real = GSL_REAL(this->data[i*this->N+j]);
                double imag = GSL_IMAG(this->data[i*this->N+j]);

                double lx, ly;

                if(i > this->N / 2)
                {
                    lx = i - this->N;
                }
                else
                {
                    lx = i;
                }
                if(j > this->N / 2)
                {
                    ly = j - this->N;
                }
                else
                {
                    ly = j;
                }

                if(i != 0 || j != 0)
                {
                    double filter = (2.0*lx*ly) / (lx*lx + ly*ly) / e_crit;
                    this->set(i,j,real*filter,imag*filter);
                }
                else
                {
                    double filter = 0.0;
                    this->set(i,j,real*filter,imag*filter);
                }
            }
        }
    }

    void reverse_filter_shear2(double e_crit)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real = GSL_REAL(this->data[i*this->N+j]);
                double imag = GSL_IMAG(this->data[i*this->N+j]);

                double lx, ly;

                if(i > this->N / 2)
                {
                    lx = i - this->N;
                }
                else
                {
                    lx = i;
                }
                if(j > this->N / 2)
                {
                    ly = j - this->N;
                }
                else
                {
                    ly = j;
                }

                if(i != 0 || j != 0)
                {
                    double filter = e_crit * (2.0*lx*ly) / (lx*lx + ly*ly);
                    //filter = filter * exp(-i*i - j*j);
                    this->set(i,j,real*filter,imag*filter);
                }
                else
                {
                    double filter = 0.0;
                    this->set(i,j,real*filter,imag*filter);
                }
            }
        }
    }

    void add_noise(double sigma, gsl_rng* rng)
    {
        for(int i = 0; i < this->N; i++)
        {
            for(int j = 0; j < this->N; j++)
            {
                double real = GSL_REAL(this->data[i*this->N+j]);
                double noise = gsl_ran_gaussian(rng, sigma);

                this->set(i,j,real+noise,0.0);

            }
        }
    }
};

#endif
