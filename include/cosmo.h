#ifndef COSMO_H
#define COSMO_H

#include <gsl/gsl_integration.h>

class cosmology
{
private:
    double Om;
    double h;
    double dh;
    gsl_integration_cquad_workspace* w;

    const static double c = 3.0e5;

    static double dc_integral(double z, void* params)
    {
        cosmology* c = (cosmology*) params;
        double Om = c->Om;

        double result = 1.0/sqrt(Om*(1.0+z)*(1.0+z)*(1.0+z)+(1.0-Om));

        return result;
    }

public:
    cosmology(double h, double Om)
    {
        this->Om = Om;
        this->h = h;

        this->dh = c / (h*100.0);

        this->w = gsl_integration_cquad_workspace_alloc(1000);
    }

    double E(double z)
    {
        double result;

        result = sqrt(Om*(1+z)*(1+z)*(1+z) + (1.0 - Om));

        return result;
    }

    double dc(double z)
    {
        double result;
        double error;
        size_t nevals;

        gsl_function f;
        f.function = &dc_integral;
        f.params = this;

        gsl_integration_cquad(&f, 0.0, z, 1e-9, 1e-7, this->w, &result, &error, &nevals);

        return dh*result;
    }

    double da(double z)
    {
        double result;

        result = dc(z) / (1.0 + z);

        return result;
    }

    ~cosmology()
    {
        if (this->w)
        {
            gsl_integration_cquad_workspace_free(this->w);
            this->w = 0;
        }
    }
};





#endif
