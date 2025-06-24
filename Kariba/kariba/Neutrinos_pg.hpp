#pragma once

#include <string>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "kariba/Radiation.hpp"

namespace kariba {

struct HetaParams {
    // eq 70 from KA08 for photons and writen as 0< x=Eg/Ep <1
    double eta;
    double eta_zero;
    double E;
    double gp_min;
    double gp_max;
    gsl_spline *spline_Jp;
    gsl_interp_accel *acc_Jp;
    std::string product;
    gsl_interp_accel *acc_ng;
    gsl_spline *spline_ng;
    double nu_min;
    double nu_max;
};

class Neutrinos_pg : public Radiation {
    //	private:
  public:
    ~Neutrinos_pg();
    Neutrinos_pg(int s1, double Emin, double Emax);

    void set_neutrinos(double gp_min, double gp_max, gsl_interp_accel *acc_Jp,
                       gsl_spline *spline_Jp, double *en_perseg,
                       double *lum_perseg, int nphot,
                       std::string outputConfiguration, std::string flavor,
                       int infosw, std::string source);
};

double Heta(double x, void *p);
double colliding_protons(gsl_spline *spline_Jp, gsl_interp_accel *acc_Jp,
                         double gp_min, double gp_max, double Ep);
double photons_jet(double eta, double Ep, gsl_spline *spline_ng,
                   gsl_interp_accel *acc_ng, double nu_min, double nu_max);
void tables_photomeson(double &s, double &delta, double &Beta,
                       std::string product, double xeta);
double PhiFunc(double eta, double eta0, double x, std::string product);

}    // namespace kariba
