#pragma once

#include <string>

#include <gsl/gsl_spline.h>

namespace kariba {

// Template class for photon/neutrino distributions

// Structures used for GSL integration
struct CyclosynEmisParams {
    double nu;
    double b;
    gsl_spline *syn;
    gsl_interp_accel *acc_syn;
    gsl_spline *eldis;
    gsl_interp_accel *acc_eldis;
};

struct CyclosynAbsParams {
    double nu;
    double b;
    gsl_spline *syn;
    gsl_interp_accel *acc_syn;
    gsl_spline *derivs;
    gsl_interp_accel *acc_derivs;
};

struct ComintParams {
    double eph;
    double ephmin;
    double ephmax;
    gsl_spline *eldis;
    gsl_interp_accel *acc_eldis;
    gsl_spline *phodis;
    gsl_interp_accel *acc_phodis;
};

struct ComfncParams {
    double game;
    double e1;
    gsl_spline *phodis;
    gsl_interp_accel *acc_phodis;
};

struct DiskObsParams {
    double tin;
    double rin;
    double nu;
};

struct DiskIcParams {
    double gamma;
    double beta;
    double tin;
    double rin;
    double rout;
    double h;
    double z;
    double nu;
};

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
    double get_volume() const { return vol; }

    double integrated_luminosity(double numin, double numax);

    void set_beaming(double theta, double speed, double doppler);
    void set_inclination(double theta);
    void set_geometry(std::string geom, double l1, double l2);
    void set_geometry(std::string geom, double l1);

    void set_counterjet(bool flag);
    void test_arrays();
};
}    // namespace kariba
