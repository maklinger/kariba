#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>

const double cee = GSL_CONST_CGSM_SPEED_OF_LIGHT;
const double emgm = GSL_CONST_CGSM_MASS_ELECTRON;
const double pmgm = GSL_CONST_CGSM_MASS_PROTON;
const double kboltz = GSL_CONST_CGSM_BOLTZMANN;
const double kboltz_kev2erg = 1.6022e-9;
const double gr_to_kev = 5.6095883571872e+29;
const double me_kev = 511.0;
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
const double erg = 6.24e11;    // 1 erg = 6.24e11 eV
const double Kpp  = 0.5;        // Inelasticity Kpp. Here is considered constant
const double barn = 1.e-24;
const double mbarn = 1.e-3 * barn;
const double sigmapp = 3.43e-26;    // pp cross section in cm2

// Template class for particle distributions
// This class contains members and methods that are used for thermal,
// non-thermal and mixed distributions

// Structures used for GSL integration
typedef struct pl_params {
    double s;
    double n;
} pl_params;

typedef struct bkn_params {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
} bkn_params;

typedef struct th_params {
    double t;
    double n;
    double m;
} th_params;

typedef struct k_params {
    double t;
    double k;
} k_params;

typedef struct injection_mixed_params {
    double s;
    double t;
    double nth;
    double npl;
    double m;
    double min;
    double max;
    double cutoff;
} injection_mixed_params;

typedef struct injection_kappa_params {
    double t;
    double k;
    double n;
    double m;
} injection_kappa_params;

typedef struct injection_pl_params {
    double s;
    double n;
    double m;
    double max;
} injection_pl_params;

typedef struct injection_bkn_params {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
    double n;
} injection_bkn_params;

class Particles {
  protected:
    int size;

    double mass_gr;    // particle mass in grams
    double
        mass_kev;    // same as above but in keV, using electrons as "reference"

    double *p;    // array of particle momenta
    double
        *ndens;    // array of number density per unit volume, per unit momentum
    double *gamma;    // array of particle kinetic energies for each momentum
    double *gdens;    // array of number density per unit volume, per unit gamma
    double *gdens_diff;    // array with differential of number density for
                           // radiation calculation

  public:
    ~Particles();

    void set_mass(double m);
    void initialize_gdens();
    void initialize_pdens();
    void gdens_differentiate();

    const double *get_p() const { return p; }
    const double *get_pdens() const { return ndens; }
    const double *get_gamma() const { return gamma; }
    const double *get_gdens() const { return gdens; }
    const double *get_gdens_diff() const { return gdens_diff; }

    double beta(int i);

    double count_particles();
    double count_particles_energy();
    double av_p();
    double av_gamma();
    double av_psq();
    double av_gammasq();

    void test_arrays();
};
#endif
