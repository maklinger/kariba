#pragma once

namespace kariba {
// Template class for particle distributions
// This class contains members and methods that are used for thermal,
// non-thermal and mixed distributions

// Structures used for GSL integration
struct PlParams {
    double s;
    double n;
};

struct BknParams {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
};

struct ThParams {
    double t;
    double n;
    double m;
};

struct KParams {
    double t;
    double k;
};

struct InjectionMixedParams {
    double s;
    double t;
    double nth;
    double npl;
    double m;
    double min;
    double max;
    double cutoff;
};

struct InjectionKappaParams {
    double t;
    double k;
    double n;
    double m;
};

struct InjectionPlParams {
    double s;
    double n;
    double m;
    double max;
};

struct InjectionBknParams {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
    double n;
};

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

    double *get_p() const { return p; }
    double *get_pdens() const { return ndens; }
    double *get_gamma() const { return gamma; }
    double *get_gdens() const { return gdens; }
    double *get_gdens_diff() const { return gdens_diff; }

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
