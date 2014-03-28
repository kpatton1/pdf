#ifndef GF_H
#define GF_H

class GF
{
public:

    virtual double g(double z) { return 1.0; }
    virtual ~GF(){}
};

class GF_LCDM : public GF
{
public:
    double g(double z) { return 1.0; }
};


#endif

