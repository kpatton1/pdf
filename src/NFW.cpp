
#include <cmath>

double NFW_projected_mass(double x)
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

