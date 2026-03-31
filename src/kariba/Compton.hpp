#pragma once

#include <gsl/gsl_spline2d.h>

#include "Radiation.hpp"

namespace kariba {

//! Class inverse Compton, inherited from Radiation.hpp
class Compton : public Radiation {
  protected:
    size_t Niter;          //!< number of IC iterations
    double tau, ypar;      //!< optical depth/comtpon Y of emitting region
    double rphot;          //!< photospheric radius when tau > 1, used to renormalize volume
    double escape_corr;    //!< escape term, used to renormalize our spectra to CompPS

    std::vector<double> log_target_energy;    //!< array of seed energies in erg
    std::vector<double> log_target_diff_spec;     //!< array of seed photon number density in log(#/erg/cm^3)
    std::vector<double>
        log_target_diff_spec_iter;    //!< array of iterated photon number density in log(#/erg/cm^3)
    std::vector<double> log_energy_iter;

    gsl_spline* seed_ph;           //!< interpolation of photon field array target_diff_spec
    gsl_interp_accel* acc_seed;    //!< accelerator for above spline

    gsl_spline* iter_ph;           //!< interpolation of photon field for multiple scatters
    gsl_interp_accel* acc_iter;    //!< accelerator of above spline

    gsl_spline2d* esc_p_sph;    //!< interpolation for escape calculation to mimic
    //!< radiative transfer
    gsl_spline2d* esc_p_cyl;    //!< interpolation for escape calculation to mimic
    //!< radiative transfer
    gsl_interp_accel* acc_tau;    //!< accelerator of above spline over tau
    gsl_interp_accel* acc_Te;     //!< accelerator of above spline over Te

    gsl_integration_workspace* w1;
    gsl_integration_workspace* w2;

  public:
    ~Compton();
    Compton(size_t size, size_t target_size);

    friend double comfnc(double logein, void* pars);
    friend double comint(double gam, void* pars);
    friend double disk_integral(double alfa, void* p);
    virtual double comintegral(size_t it, double blim, double ulim, double nu, double numin, double numax,
                       gsl_spline* eldis, gsl_interp_accel* acc_eldis);
    virtual void compton_spectrum(double gmin, double gmax, gsl_spline* eldis, gsl_interp_accel* acc_eldis);

    virtual void set_target_energy_array(const std::vector<double>& new_target_energy);
    virtual void set_target_frequency_array(const std::vector<double>& new_target_frequency);

    virtual void add_target_diff_spec(const std::vector<double>& new_target_diff_spec);
    virtual void add_target_energy_density(
      const std::vector<double>& new_target_energy,
      const std::vector<double>& new_target_energy_density
    );
    virtual void add_target_number_density(
      const std::vector<double>& new_target_energy,
      const std::vector<double>& new_target_number_density
    );
    virtual void add_target_number_density_on_freq_grid(
      const std::vector<double>& new_target_frequency,
      const std::vector<double>& new_target_number_density
    );

    
    virtual void cyclosyn_seed(const std::vector<double>& syn_frequencies, const std::vector<double>& syn_number_densities);
    virtual void bb_seed_k(double Urad, double Tbb);
    virtual void bb_seed_kev(double Urad, double Tbb);
    virtual void shsdisk_seed(double tin, double rin, double rout,
                      double h, double z);

    virtual void set_niter(double nu0, double Te);
    virtual void set_niter(size_t n);
    virtual void set_tau(double n, double gam);
    virtual void set_tau(double _tau);
    virtual void set_frequency(double numin, double numax);
    virtual void set_escape(double escape);

    virtual std::vector<double> get_target_energy();
    virtual std::vector<double> get_target_diff_spec();

    virtual double get_tau() const { return tau; };

    virtual double get_ypar() const { return ypar; };

    friend double comfnc(double ein, void* p);
    friend double comint(double gam, void* p);
    friend double disk_integral(double alfa, void* p);
};

}    // namespace kariba
