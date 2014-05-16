
#include "cosmo.h"

#include <cmath>
#include <iostream>

const double G = 4.302e-9; // G in units of Mpc Msolar^-1 (km/s)^2
const double c = 3.0e5; // c in unites of km/s

int main()
{
    double h = 0.7;
    double Om = 0.30;
    double pc = 3.0 * 100.0 * 100.0 * h * h / (8.0 * M_PI * G);
    double p = Om * pc;

    cosmology* cosmo = new cosmology(h, Om);

    double z_lens = 0.2;

    double z_source = 0.6;

    double da_lens = cosmo->da(z_lens);
    double da_source = cosmo->da(z_source);

    std::cout << "da_lens: " << da_lens << "  da_source: " << da_source << std::endl;

    double n_gal = 10.0; // galaxies per sq arcmin

    std::cout << "n_gal: " << n_gal << " per sq arcmin" << std::endl;

    n_gal = n_gal * 60.0 * 60.0 * 180.0 * 180.0 / M_PI / M_PI; // galaxies per sq radian

    std::cout << "n_gal: " << n_gal << " per sq radian" << std::endl;

    n_gal = n_gal / da_lens / da_lens; // galaxies per sq Mpc

    std::cout << "n_gal: " << n_gal << " per sq Mpc" << std::endl;

    double e_crit = c*c*da_source/(4.0 * M_PI * G * (da_source - da_lens)*da_lens);

    std::cout << "e_crit/p: " << e_crit/p << std::endl;

    return 0;
}
