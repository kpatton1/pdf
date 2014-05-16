#include <cmath>

#include <iostream>
#include <fstream>

#include "gsl/gsl_interp.h"

int main()
{
    double x_array[100000];
    double y_array[100000];
    double xy_array[100000];
    double x2y_array[100000];

    int i = 0;

    while (std::cin.good() && i < 100000)
    {
        double x;
        double y;
        std::cin >> x >> y;

        if (!std::cin.good())
        {
            break;
        }

        x = exp(x);
        y = y / x;

        if (i > 0 && x <= x_array[i - 1])
        {
            continue;
        }

        x_array[i] = x;
        y_array[i] = y;
        xy_array[i] = x * y;
        x2y_array[i] = x * x * y;

        i++;

    }

    //std::cout << i << std::endl;

    gsl_interp* interp_y = gsl_interp_alloc(gsl_interp_cspline, i);
    gsl_interp* interp_xy = gsl_interp_alloc(gsl_interp_cspline, i);
    gsl_interp* interp_x2y = gsl_interp_alloc(gsl_interp_cspline, i);

    double x_min = x_array[0];
    double x_max = x_array[i - 1];

    gsl_interp_init(interp_y, x_array, y_array, i);
    gsl_interp_init(interp_xy, x_array, xy_array, i);
    gsl_interp_init(interp_x2y, x_array, x2y_array, i);

    gsl_interp_accel* accel = gsl_interp_accel_alloc();

    //std::cout << integral_value << std::endl;

    for (int j = 0; j < i; j++)
    {
        double integral_y = gsl_interp_eval_integ(interp_y, x_array, y_array, x_min, x_array[j], accel);
        double integral_xy = gsl_interp_eval_integ(interp_xy, x_array, xy_array, x_min, x_array[j], accel);
        double integral_x2y = gsl_interp_eval_integ(interp_x2y, x_array, x2y_array, x_min, x_array[j], accel);

        std::cout << x_array[j] << " " << y_array[j] << " " << integral_y << " " << integral_xy << " " << integral_x2y << std::endl;
    }

    double integral_y = gsl_interp_eval_integ(interp_y, x_array, y_array, x_min, x_max, accel);
    double integral_xy = gsl_interp_eval_integ(interp_xy, x_array, xy_array, x_min, x_max, accel);
    double integral_x2y = gsl_interp_eval_integ(interp_x2y, x_array, x2y_array, x_min, x_max, accel);

    std::cerr << "integral y(x): " << integral_y << "  integral y(x)*x: " << integral_xy << "  integral y(x)*x*x: " << integral_x2y << std::endl;
    std::cerr << "mean: " << integral_xy << "  std: " << integral_x2y - integral_xy*integral_xy << std::endl;

    gsl_interp_accel_free(accel);
    gsl_interp_free(interp_y);
    gsl_interp_free(interp_xy);
    gsl_interp_free(interp_x2y);

    return 0;
}
