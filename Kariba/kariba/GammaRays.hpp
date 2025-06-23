/*************************************************************************************************************
Gamma-rays from neutral pion decay, products of inelastic pp and pγ collisions
*************************************************************************************************************/
#ifndef GAMMARAYS_HPP
#define GAMMARAYS_HPP

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "kariba/Radiation.hpp"

namespace kariba {

#ifndef PHOTOMESON_TABLES
#define PHOTOMESON_TABLES

//----------------------- γ rays -----------------------//
static const double etagTable[22] = {1.1, 1.2,  1.3,  1.4, 1.5,  1.6,  1.7, 1.8,
                                     1.9, 2.0,  3.0,  4.0, 5.0,  6.0,  7.0, 8.0,
                                     9.0, 10.0, 20.0, 30., 40.0, 100.0};
static const double sgTable[22] = {
    0.0768, 0.106,  0.182, 0.201,  0.219,  0.216, 0.233, 0.233,
    0.248,  0.244,  0.188, 0.131,  0.120,  0.107, 0.102, 0.0932,
    0.0838, 0.0761, 0.107, 0.0928, 0.0722, 0.0479};
static const double deltagTable[22] = {
    0.544, 0.540, 0.750, 0.791, 0.788, 0.831, 0.839, 0.825, 0.805, 0.779, 1.23,
    1.82,  2.05,  2.19,  2.23,  2.29,  2.37,  2.43,  2.27,  2.33,  2.42,  2.59};
static const double BetagTable[22] = {
    2.86e-19, 2.24e-18, 5.61e-18, 1.02e-17, 1.60e-17, 2.23e-17,
    3.10e-17, 4.07e-17, 5.30e-17, 6.74e-17, 1.51e-16, 1.24e-16,
    1.37e-16, 1.62e-16, 1.71e-16, 1.78e-16, 1.84e-16, 1.93e-16,
    4.74e-16, 7.70e-16, 1.06e-15, 2.73e-15};

#endif

typedef struct Hetag_params {
    // eq 70 from KA08 for photons and writen as 0< x=Eg/Ep <1
    double eta;
    double eta_zero;
    double Eg;
    double gp_min;
    double gp_max;
    gsl_spline *spline_Jp;
    gsl_interp_accel *acc_Jp;
    gsl_interp_accel *acc_ng;
    gsl_spline *spline_ng;
    double nu_min;
    double nu_max;
} Hetag_params;

class Grays : public Radiation {
    //	private:
  public:
    ~Grays();
    Grays(int s1, double numin, double numax);

    // Method to set the gamma-rays from pp inelastic interactions. p: pspec_p,
    // ntot_prot: total proton number density of the jet segment,
    // ntargets: the number density of external proton density (companion etc),
    // plfrac: plfrac_p
    void set_grays_pp(double p, double gammap_min, double gammap_max,
                      double ntot_prot, double ntargets, double plfrac,
                      gsl_interp_accel *acc_Jp, gsl_spline *spline_Jp);

    void set_grays_pg(double gp_min, double gp_max, gsl_interp_accel *acc_Jp,
                      gsl_spline *spline_Jp, double *nu_per_seg,
                      double *ng_per_seg, int ne);
};

// Adds up in the lum_perseg the target photon luminosity (in erg/sec/Hz)
void sum_photons(int nphot, double *en_perseg, double *lum_perseg, int ntarg,
                 const double *targ_en, const double *targ_lum);
void sum_photons(int nphot, const double *en_perseg, double *lum_perseg,
                 int ntarg, const double *targ_en, const double *targ_lum);

// funtions for γ rays from pγ
double Hetag(double x, void *p);

// The following are common for γ rays/electrons/neutrinos from pp:
double set_ntilde(double p);
double target_protons(double ntot_prot, double ntargets, double plfrac);
double sigma_pp(double Ep);
double proton_dist(double gpmin, double Ep, double Epcode_max,
                   gsl_spline *spline_Jp, gsl_interp_accel *acc_Jp);
double gspec_pp(double Ep, double y);

// The following are common for γ rays/electrons/neutrinos from pγ:
double colliding_protons(gsl_spline *spline_Jp, gsl_interp_accel *acc_Jp,
                         double gp_min, double gp_max, double Ep);
double photons_jet(double eta, double Ep, gsl_spline *spline_ng,
                   gsl_interp_accel *acc_ng, double nu_min, double nu_max);
void tables_photomeson_gamma(double &s, double &delta, double &Beta,
                             double xeta);
double PhiFunc_gamma(double eta, double eta0, double x);

}    // namespace kariba

#endif
