#ifndef RADIATION_HPP
#define RADIATION_HPP

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <string>

const double kpc = 1e3 * GSL_CONST_CGSM_PARSEC;
const double cee = GSL_CONST_CGSM_SPEED_OF_LIGHT;
const double emgm = GSL_CONST_CGSM_MASS_ELECTRON;
const double pmgm = GSL_CONST_CGSM_MASS_PROTON;
const double kboltz = GSL_CONST_CGSM_BOLTZMANN;
const double kboltz_kev2erg = 1.6022e-9;    // Boltzman constant in keV/erg
const double emerg = GSL_CONST_CGSM_MASS_ELECTRON * pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0);
const double pi = M_PI;
const double charg = 4.8e-10;
const double sigtom = GSL_CONST_CGSM_THOMSON_CROSS_SECTION;
const double herg = GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
const double hkev = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 6.2415e8;
const double mjy = 1.e-26;
const double re0 = 2.81794e-13;
const double gconst = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT;
const double sbconst = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT;
const double aconst = 7.56e-15;
const double msun = GSL_CONST_CGSM_SOLAR_MASS;
const double erg = 6.24e11;               // 1 erg = 6.24e11 eV
const double mprotTeV = 938.272046e-6;    // mass of proton in TeV/c^2
const double mpionTeV = 139.57e-6;        // mass of pion in TeV/c^2
const double Kpp = 0.5;     // Inelasticity Kpp. Here is considered constant.
const double Kpi = 0.17;    // fraction of E_kinetic of proton transferred to neutrinos
const double hbar = herg / (2.0 * pi);    // h bar
const double barn = 1.0e-24;
const double mbarn = 1.e-3 * barn;

// Template class for photon/neutrino distributions

// Structures used for GSL integration
typedef struct cyclosyn_emis_params {
    double nu;
    double b;
    gsl_spline *syn;
    gsl_interp_accel *acc_syn;
    gsl_spline *eldis;
    gsl_interp_accel *acc_eldis;
} cyclosyn_emis_params;

typedef struct cyclosyn_abs_params {
    double nu;
    double b;
    gsl_spline *syn;
    gsl_interp_accel *acc_syn;
    gsl_spline *derivs;
    gsl_interp_accel *acc_derivs;
} cyclosyn_abs_params;

typedef struct comint_params {
    double eph;
    double ephmin;
    double ephmax;
    gsl_spline *eldis;
    gsl_interp_accel *acc_eldis;
    gsl_spline *phodis;
    gsl_interp_accel *acc_phodis;
} comint_params;

typedef struct comfnc_params {
    double game;
    double e1;
    gsl_spline *phodis;
    gsl_interp_accel *acc_phodis;
} comfnc_params;

typedef struct disk_obs_params {
    double tin;
    double rin;
    double nu;
} disk_obs_params;

typedef struct disk_ic_params {
    double gamma;
    double beta;
    double tin;
    double rin;
    double rout;
    double h;
    double z;
    double nu;
} disk_ic_params;

class Radiation {
  protected:
    int size;                // size of arrays
    double *en_phot;         // array of photon energies
    double *num_phot;        // array of number of photons in units of erg/s/Hz
    double *en_phot_obs;     // same as above but in observer frame
    double *num_phot_obs;    // same as above but in observer frame

    double r, z;             // Dimensions of emitting region
    double vol;              // Volume of emitting region
    double beta;             // speed of the emitting region
    double dopfac, angle;    // Viewing angle/Doppler factor of emitting region
    double dopnum;           // Doppler boosting exponent, depends on geometry
    bool counterjet;    // boolean switch if user wants to include counterjet
                        // emission
    std::string geometry;    // string to track geometry of emitting region

  public:
    double *get_energy() const { return en_phot; }
    double *get_nphot() const { return num_phot; }
    double *get_energy_obs() const { return en_phot_obs; }
    double *get_nphot_obs() const { return num_phot_obs; }
    int get_size() const { return size; }
    double get_volume() const         { return vol; }

    double integrated_luminosity(double numin, double numax);

    void set_beaming(double theta, double speed, double doppler);
    void set_inclination(double theta);
    void set_geometry(std::string geom, double l1, double l2);
    void set_geometry(std::string geom, double l1);

    void set_counterjet(bool flag);
    void test_arrays();
};

#endif
