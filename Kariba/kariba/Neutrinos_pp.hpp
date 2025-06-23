#ifndef NEUTRINOS_PP_HPP
#define NEUTRINOS_PP_HPP

#include <cstring>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "kariba/Radiation.hpp"

namespace kariba {

class Neutrinos_pp : public Radiation {
    //	private:
  public:
    ~Neutrinos_pp();
    Neutrinos_pp(int s1, double Emin, double Emax);

    void set_neutrinos_pp(double p, double gammap_min, double gammap_max,
                          double ntot_prot, double nwind, double plfrac,
                          gsl_interp_accel *acc_Jp, gsl_spline *spline_Jp,
                          std::string outputConfiguration, std::string flavor,
                          int infosw, std::string source);
};

double multiplicity(double pspec);    // in Electrons.cpp
double target_protons(double ntot_prot, double nwind,
                      double plfrac);    // in Gamma_rays.cpp
double prob();                           // in Electrons.cpp
double sigma_pp(double Ep);              // in Gamma_rays.cpp
double proton_dist(double gpmin, double Ep, double Epcode_max,
                   gsl_spline *spline_Jp, gsl_interp_accel *acc_Jp);
// in Gamma_rays.cpp

double distr_pp(double lEv, double lEpi, std::string flavor);
double secondary_spectrum(double Ep, double y, std::string flavor);

double prob_fve();

}    // namespace kariba

#endif
