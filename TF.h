#ifndef TF_H
#define TF_H

#include <cmath>

class TF
{
public:
    virtual double t(double k) { return 0.0; };
    virtual double t_k2(double k) { return 0.0; };

};

class TF_BBKS : public TF
{
private:
    double Om;
    double h;
    
    double factor;

public:
    TF_BBKS(double Om, double h)
    {
        this->Om = Om;
        this->h = h;
        this->factor = Om * h * h;
    };

    double t(double k)
    {
        double q = k / factor;
        
        double result = log(1.0 + 2.34*q) * pow(1.0 + 3.89*q + 16.1*16.1*q*q + 5.46*5.46*5.46*q*q*q + 6.71*6.71*6.71*6.71*q*q*q*q, -0.25) / (2.34*q);
        
        return result;
    };
    
    double t_k2(double k)
    {
        double q = k / factor;
        
        double result = factor * factor * log(1.0 + 2.34*q) * q * pow(1.0 + 3.89*q + 16.1*16.1*q*q + 5.46*5.46*5.46*q*q*q + 6.71*6.71*6.71*6.71*q*q*q*q, -0.25) / (2.34);
       
        return result;
    };
};

#endif

