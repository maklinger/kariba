#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "kariba/Particles.hpp"
#include "kariba/constants.hpp"

namespace kariba {

Particles::Particles(size_t size)
    : p(size, 0.0), ndens(size, 0.0), gamma(size, 0.0), gdens(size, 0.0), gdens_diff(size, 0.0) {}

//! Simple numerical integrals /w trapeze method
double Particles::count_particles() {
    double temp = 0.0;
    for (size_t i = 0; i < ndens.size() - 1; i++) {
        temp = temp + (1. / 2.) * (p[i + 1] - p[i]) * (ndens[i + 1] + ndens[i]);
    }
    return temp;
}

double Particles::count_particles_energy() {
    double temp = 0.0;
    for (size_t i = 0; i < gamma.size() - 1; i++) {
        temp = temp + (1. / 2.) * (gamma[i + 1] - gamma[i]) * (gdens[i + 1] + gdens[i]);
    }
    return temp;
}

double Particles::av_p() {
    double temp = 0.;
    for (size_t i = 0; i < p.size() - 1; i++) {
        temp = temp + (1. / 2.) * (p[i + 1] - p[i]) * (p[i + 1] * ndens[i + 1] + p[i] * ndens[i]);
    }
    return temp / count_particles();
}

double Particles::av_gamma() {
    return pow(pow(av_p() / (mass_gr * constants::cee), 2.) + 1., 1. / 2.);
}

double Particles::av_psq() {
    double temp = 0.;
    for (size_t i = 0; i < p.size() - 1; i++) {
        temp = temp + (1. / 2.) * (p[i + 1] - p[i]) *
                          (pow(p[i + 1], 2.) * ndens[i + 1] + pow(p[i], 2.) * ndens[i]);
    }
    return temp / count_particles();
}

double Particles::av_gammasq() {
    return pow(av_psq() / pow(mass_gr * constants::cee, 2.) + 1., 1. / 2.);
}

//! Methods to set up energy space number density, as a function of momentum
//! space number density
void Particles::initialize_gdens() {
    for (size_t i = 0; i < gdens.size(); i++) {
        gdens[i] =
            ndens[i] * gamma[i] * mass_gr * constants::cee / (pow(pow(gamma[i], 2.) - 1., 1. / 2.));
    }
}

//! Same as above but the other way around
void Particles::initialize_pdens() {
    for (size_t i = 0; i < gdens.size(); i++) {
        ndens[i] = gdens[i] * p[i] /
                   (pow(mass_gr * constants::cee, 2.) *
                    pow(pow(p[i] / (mass_gr * constants::cee), 2.) + 1., 1. / 2.));
    }
}

void Particles::gdens_differentiate() {
    std::vector<double> temp;
    size_t size = gdens.size();

    for (size_t i = 0; i < size; i++) {
        temp.push_back(gdens[i] / (pow(gamma[i], 1.)));
    }

    for (size_t i = 0; i < size - 1; i++) {
        gdens_diff[i] = (temp[i + 1] - temp[i]) /
                        (mass_gr * pow(constants::cee, 2.) * (gamma[i + 1] - gamma[i]));
    }

    gdens_diff[size - 1] = gdens_diff[size - 2];
}

void Particles::set_mass(double m) {
    mass_gr = m;
    mass_kev = m * constants::gr_to_kev;
}

//! simple method to check arrays; only meant for debugging
void Particles::test_arrays() {
    for (size_t i = 0; i < p.size(); i++) {
        std::cout << p[i] << " " << gamma[i] << " " << ndens[i] << " " << ndens[i] * p[i]
                  << std::endl;
    }
}

}    // namespace kariba
