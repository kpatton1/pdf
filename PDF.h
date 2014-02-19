#ifndef PDF_H
#define PDF_H

#include "gsl/gsl_integration.h"

class PDF
{
private:

public:
    

};

double gsl_area_integral(double lnm, void* params);
double gsl_pdf_integral(double lnm, void* params);

struct pdf_var_struct
{
    void* pdf;
    double z;
};

struct pdf_var_struct2
{
    void* pdf;
    double lne;
    double z;
};

class PDF_HALO : public PDF
{
private:
    MF* mf;
    Halo* halo;

    //gsl_integration_workspace* w;
    
    gsl_integration_cquad_workspace* w;

public:
    friend double gsl_area_integral(double lnm, void* params);
    friend double gsl_pdf_integral(double lnm, void* params);
    
    PDF_HALO(MF* mf, Halo* halo)
    {
        this->mf = mf;
        this->halo = halo;
        
        //this->w = gsl_integration_workspace_alloc(1000000);
        this->w = gsl_integration_cquad_workspace_alloc(10000);
    };
    
    double covering(double lnm_min, double lnm_max, double z)
    {
        double result;
        double error;
        size_t nevals;
        
        gsl_function f;
        f.function = &gsl_area_integral;
        
        pdf_var_struct params;
        params.pdf = this;
        params.z = z;
        
        f.params = &params;
        
        //gsl_integration_qag(&f, lnm_min, lnm_max, 1e-12, 1e-8, 1000000, GSL_INTEG_GAUSS51, w, &result, &error);
        gsl_integration_cquad(&f, lnm_min, lnm_max, 1e-9, 1e-7, this->w, &result, &error, &nevals);
        
        return (result * M_PI);
    };
    
    double pdf(double lne, double lnm_min, double lnm_max, double z)
    {
        double result;
        double error;
        size_t nevals;
        
        gsl_function f;
        f.function = &gsl_pdf_integral;
        
        pdf_var_struct2 params;
        params.pdf = this;
        params.lne = lne;
        params.z = z;
        
        f.params = &params;
        
        //gsl_integration_qag(&f, lnm_min, lnm_max, 1e-12, 1e-8, 1000000, GSL_INTEG_GAUSS51, w, &result, &error);
        gsl_integration_cquad(&f, lnm_min, lnm_max, 1e-9, 1e-7, this->w, &result, &error, &nevals);
        
        return (result * M_PI);
    };
    
    ~PDF_HALO()
    {
        if(this->w)
        {
            //gsl_integration_workspace_free(this->w);
            gsl_integration_cquad_workspace_free(this->w);
        }
    }
};

double gsl_area_integral(double lnm, void* params)
{
    pdf_var_struct* p = (pdf_var_struct*) params;
    PDF_HALO* pdf = (PDF_HALO*) p->pdf;
    double z = p->z;
    double dn_dlnm = pdf->mf->dn_dlnm(lnm, z);

    //double r200 = pow(m_p * 3.0 / (800.0 * M_PI), 1.0/3.0);
    double r200 = exp(1.0/3.0 * (lnm + log(3.0 / (800.0 * M_PI))));

    double result = dn_dlnm * r200 * r200;
    
    return result;
}

double gsl_pdf_integral(double lnm, void* params)
{
    pdf_var_struct2* p = (pdf_var_struct2*) params;
    PDF_HALO* pdf = (PDF_HALO*) p->pdf;
    double lne = p->lne;
    double pdf_halo = pdf->halo->PDF(lne, lnm, 10.0);
    double z = p->z;
    double dn_dlnm = pdf->mf->dn_dlnm(lnm, z);

    //double r200 = pow(m_p * 3.0 / (800.0 * M_PI), 1.0/3.0);
    double r200 = exp(1.0/3.0 * (lnm + log(3.0 / (800.0 * M_PI))));

    double result = pdf_halo * dn_dlnm * r200 * r200;
    
    return result;
}

#endif


