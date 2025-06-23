#ifndef ELECTRONS_HPP
#define ELECTRONS_HPP

#include <gsl/gsl_spline.h>

namespace kariba {

// functions for electrons from pp
double multiplicity(double pspec);
double prob();
double elec_dist_pp(double zen, double w);
double elec_spec_pp(double Ep, double y);

double target_protons(double ntot_prot, double nwind, double plfrac);
double proton_dist(double gpmin, double Ep, double Epcode_max,
                   gsl_spline *spline_Jp, gsl_interp_accel *acc_Jp);

// function for electrons from γγ annihilation
double production_rate(double ge, double x);
}    // namespace kariba

#endif
