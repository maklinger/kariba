#pragma once

#include "Radiation.hpp"

namespace kariba {

//! Class synchrotron photons, inherited from Radiation.hpp
class Cyclosyn : public Radiation {
  protected:
    double bfield;     // Magnetic field in emitting region
    double mass_gr;    // Mass of the emitting particle
    gsl_spline* syn_f;
    gsl_interp_accel* syn_acc;
    std::vector<double> cyclosyn_absorption_rate;

  public:
    ~Cyclosyn();
    Cyclosyn(size_t size);

    virtual double emis_integral(double nu, double gmin, double gmax, gsl_spline* eldis,
                                 gsl_interp_accel* acc_eldis);
    virtual double abs_integral(double nu, double gmin, double gmax, gsl_spline* eldis_diff,
                                gsl_interp_accel* acc_eldis_diff);

    virtual void cycsyn_spectrum(double gmin, double gmax, gsl_spline* eldis,
                                 gsl_interp_accel* acc_eldis, gsl_spline* eldis_diff,
                                 gsl_interp_accel* acc_eldis_diff);

    virtual double nu_syn(double gamma);
    virtual double nu_syn();

    virtual void set_frequency(double numin, double numax);
    virtual void set_bfield(double b);
    virtual void set_mass(double mass);

    virtual void test();

    const std::vector<double>& get_cyclosyn_absorption_rate() const { return cyclosyn_absorption_rate; }

    friend double cyclosyn_emis(double gamma, void* pars);
    friend double cyclosyn_abs(double gamma, void* pars);
};

}    // namespace kariba
