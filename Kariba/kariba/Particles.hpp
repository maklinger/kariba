#pragma once

namespace kariba {
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

}    // namespace kariba
