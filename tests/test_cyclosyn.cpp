#include "doctest.h"

#include <cmath>
#include <gsl/gsl_spline.h>
#include <kariba/Compton.hpp>
#include <kariba/Cyclosyn.hpp>
#include <kariba/Powerlaw.hpp>
#include <kariba/ShSDisk.hpp>
#include <kariba/Thermal.hpp>
#include <kariba/constants.hpp>

namespace karcst = kariba::constants;

TEST_CASE("Integration tests - Complete workflows") {
    SUBCASE("Single zone jet model workflow") {
        // This test replicates the singlezone example workflow

        size_t nel = 50;    // Smaller arrays for faster testing
        size_t nfreq = 50;

        // Model parameters (simplified from example)
        double Mbh = 6.5e9;
        double B = 1.5e-3;
        double Rg = karcst::gconst * Mbh * karcst::msun / karcst::cee_cee;
        double R = 626.0 * Rg;
        double n = 9.5e-3;
        double gmin = 4.1e3;
        double pmin = std::sqrt(gmin * gmin - 1.0) * karcst::emgm * karcst::cee;
        double gmax = 6.4e4;    // Reduced for faster testing
        double p = 3.03;

        // Set up electron distribution
        kariba::Powerlaw electrons(nel);
        electrons.set_p(pmin, gmax);
        electrons.set_pspec(p);
        electrons.set_norm(n);
        electrons.set_ndens();

        CHECK(!electrons.get_p().empty());
        CHECK(electrons.count_particles() > 0.0);

        // Set up GSL interpolation
        gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
        gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen, nel);
        gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
        gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen, nel);

        gsl_spline_init(spline_eldis, electrons.get_gamma().data(), electrons.get_gdens().data(),
                        nel);
        gsl_spline_init(spline_deriv, electrons.get_gamma().data(),
                        electrons.get_gdens_diff().data(), nel);

        // Calculate synchrotron emission
        kariba::Cyclosyn syncro(nfreq);
        syncro.set_frequency(1e10, 1e18);    // Reduced range for testing
        syncro.set_bfield(B);
        syncro.set_beaming(0.0, 0.0, 1.0);    // No beaming for simplicity
        syncro.set_geometry("sphere", R);

        double gmin_actual = electrons.get_gamma()[0];
        double gmax_actual = electrons.get_gamma()[nel - 1];

        CHECK_NOTHROW(syncro.cycsyn_spectrum(gmin_actual, gmax_actual, spline_eldis, acc_eldis,
                                             spline_deriv, acc_deriv));

        // Verify synchrotron spectrum
        auto syn_energy = syncro.get_energy();
        auto syn_flux = syncro.get_nphot();

        bool has_syn_emission = false;
        for (size_t i = 0; i < nfreq; i++) {
            CHECK(syn_energy[i] > 0.0);
            if (syn_flux[i] > 0.0) {
                has_syn_emission = true;
            }
        }
        CHECK(has_syn_emission);

        // Calculate SSC emission
        kariba::Compton ssc(nfreq, nfreq);
        ssc.set_frequency(1e18, 1e24);    // Higher energy range
        ssc.set_beaming(0.0, 0.0, 1.0);
        ssc.set_geometry("sphere", R);
        ssc.set_tau(n, electrons.av_gamma() * 511.0);

        // Use synchrotron photons as seed
        ssc.cyclosyn_seed(syncro.get_energy(), syncro.get_nphot());

        CHECK_NOTHROW(ssc.compton_spectrum(gmin_actual, gmax_actual, spline_eldis, acc_eldis));

        // Verify SSC spectrum
        auto ssc_energy = ssc.get_energy();
        auto ssc_flux = ssc.get_nphot();

        bool has_ssc_emission = false;
        for (size_t i = 0; i < nfreq; i++) {
            CHECK(ssc_energy[i] >= syn_energy[nfreq - 1]);    // SSC should be higher energy
            if (ssc_flux[i] > 0.0) {
                has_ssc_emission = true;
            }
        }
        CHECK(has_ssc_emission);

        // Calculate physical quantities
        double Ue = electrons.av_gamma() * n * karcst::emgm * karcst::cee_cee;
        double Ub = B * B / (8.0 * karcst::pi);
        double equip = Ue / Ub;

        CHECK(Ue > 0.0);
        CHECK(Ub > 0.0);
        CHECK(equip > 0.0);

        // Cleanup GSL objects
        gsl_spline_free(spline_eldis);
        gsl_interp_accel_free(acc_eldis);
        gsl_spline_free(spline_deriv);
        gsl_interp_accel_free(acc_deriv);
    }
}
