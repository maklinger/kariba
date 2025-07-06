#pragma once

#include <vector>

namespace kariba {

//! Structure used for GSL integration
struct PlParams {
    double s;
    double n;
};

//! Structure used for GSL integration
struct BknParams {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
};

//! Structure used for GSL integration
struct ThParams {
    double t;
    double n;
    double m;
};

//! Structure used for GSL integration
struct KParams {
    double t;
    double k;
};

//! Structure used for GSL integration
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

//! Structure used for GSL integration
struct InjectionKappaParams {
    double t;
    double k;
    double n;
    double m;
};

//! Structure used for GSL integration
struct InjectionPlParams {
    double s;
    double n;
    double m;
    double max;
};

//! Structure used for GSL integration
struct InjectionBknParams {
    double s1;
    double s2;
    double brk;
    double max;
    double m;
    double n;
};

//! Template class for particle distributions
//! This class contains members and methods that are used for thermal,
//! non-thermal and mixed distributions
class Particles {
  protected:
    size_t size;

    double mass_gr;     //!< particle mass in grams
    double mass_kev;    //!< same as above but in keV, using electrons as "reference"

    std::vector<double> p;        //!< array of particle momenta
    std::vector<double> ndens;    //!< array of number density per unit volume, per unit momentum
    std::vector<double> gamma;    //!< array of particle kinetic energies for each momentum
    std::vector<double> gdens;    //!< array of number density per unit volume, per unit gamma
    std::vector<double> gdens_diff;    //!< array with differential of number
                                       //!< density for radiation calculation

  public:
    Particles(size_t size);

    void set_mass(double m);
    void initialize_gdens();
    void initialize_pdens();
    void gdens_differentiate();

    const std::vector<double> &get_p() const { return p; }

    const std::vector<double> &get_pdens() const { return ndens; }

    const std::vector<double> &get_gamma() const { return gamma; }

    const std::vector<double> &get_gdens() const { return gdens; }

    const std::vector<double> &get_gdens_diff() const { return gdens_diff; }

    // std::vector<double>& get_p() { return p; }
    // std::vector<double>& get_pdens() { return ndens; }
    // std::vector<double>& get_gamma() { return gamma; }
    // std::vector<double>& get_gdens() { return gdens; }
    // std::vector<double>& get_gdens_diff() { return gdens_diff; }

    // double beta(int i);  // todo: not implemented!

    double count_particles();
    double count_particles_energy();
    double av_p();
    double av_gamma();
    double av_psq();
    double av_gammasq();

    void test_arrays();
};

}    // namespace kariba
