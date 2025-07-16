#include "doctest.h"

#include <cmath>
#include <gsl/gsl_spline.h>
#include <kariba/Bknpower.hpp>
#include <kariba/Compton.hpp>
#include <kariba/Cyclosyn.hpp>
#include <kariba/Kappa.hpp>
#include <kariba/Mixed.hpp>
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
        gsl_interp_accel* acc_eldis = gsl_interp_accel_alloc();
        gsl_spline* spline_eldis = gsl_spline_alloc(gsl_interp_steffen, nel);
        gsl_interp_accel* acc_deriv = gsl_interp_accel_alloc();
        gsl_spline* spline_deriv = gsl_spline_alloc(gsl_interp_steffen, nel);

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
        const std::vector<double>& syn_energy = syncro.get_energy();
        const std::vector<double>& syn_flux = syncro.get_nphot();

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
        const std::vector<double>& ssc_energy = ssc.get_energy();
        const std::vector<double>& ssc_flux = ssc.get_nphot();

        bool has_ssc_emission = false;
        for (size_t i = 0; i < nfreq; i++) {
            CHECK(ssc_energy[i] >= syn_energy[nfreq - 1]);    // SSC should be higher/equal energy
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

    SUBCASE("Corona Comptonization workflow") {
        // This test replicates the corona example workflow

        size_t nel = 50;
        size_t nfreq = 50;

        // Parameters
        double Mbh = 10.0;
        double Rg = karcst::gconst * Mbh * karcst::msun / karcst::cee_cee;
        double Rin = 10.0 * Rg;
        double Rout = 1e4 * Rg;
        double Ldisk = 1e-4;

        // Set up disk
        kariba::ShSDisk disk;
        disk.set_mbh(Mbh);
        disk.set_rin(Rin);
        disk.set_rout(Rout);
        disk.set_luminosity(Ldisk);
        disk.set_inclination(0.0);

        CHECK_NOTHROW(disk.disk_spectrum());

        // Verify disk spectrum
        CHECK(!disk.get_energy().empty());
        CHECK(!disk.get_nphot().empty());

        double disk_lum = disk.integrated_luminosity(1e12, 1e20);
        CHECK(disk_lum > 0.0);

        // Set up corona parameters
        double Te = 90.0;    // keV
        double Tau = 0.76;
        double R_corona = 75.0 * Rg;
        double ndens = Tau / (karcst::sigtom * R_corona);

        // Set up thermal electrons
        kariba::Thermal electrons(nel);
        electrons.set_temp_kev(Te);
        electrons.set_p();
        electrons.set_norm(ndens);
        electrons.set_ndens();

        CHECK(electrons.count_particles() > 0.0);

        // Set up GSL interpolation
        gsl_interp_accel* acc_eldis = gsl_interp_accel_alloc();
        gsl_spline* spline_eldis = gsl_spline_alloc(gsl_interp_steffen, nel);

        gsl_spline_init(spline_eldis, electrons.get_gamma().data(), electrons.get_gdens().data(),
                        nel);

        // Set up Comptonization
        kariba::Compton ic(nfreq, 50);    // 50 is disk spectrum size
        ic.set_frequency(1e16, 1e20);     // Reduced range for testing
        ic.set_beaming(0.0, 0.0, 1.0);
        ic.set_geometry("sphere", R_corona);
        ic.set_tau(ndens, Te);
        ic.set_niter(5);    // Reduced iterations for testing

        // Use disk as seed photon field
        ic.shsdisk_seed(disk.get_energy(), disk.tin(), Rin, Rout, disk.hdisk(), 0.0);

        double gmin = electrons.get_gamma()[0];
        double gmax = electrons.get_gamma()[nel - 1];

        CHECK_NOTHROW(ic.compton_spectrum(gmin, gmax, spline_eldis, acc_eldis));

        // Verify Compton spectrum
        const std::vector<double>& ic_energy = ic.get_energy();
        const std::vector<double>& ic_flux = ic.get_nphot();

        bool has_ic_emission = false;
        for (size_t i = 0; i < nfreq; i++) {
            CHECK(ic_energy[i] > 0.0);
            if (ic_flux[i] > 0.0) {
                has_ic_emission = true;
            }
        }
        CHECK(has_ic_emission);

        // Check that Compton spectrum extends to higher energies than disk
        const std::vector<double>& disk_energy = disk.get_energy();
        CHECK(ic_energy[nfreq - 1] > disk_energy[49]);    // IC should extend higher

        // Cleanup
        gsl_spline_free(spline_eldis);
        gsl_interp_accel_free(acc_eldis);
    }

    SUBCASE("Particle distribution comparison workflow") {
        // This test replicates the particles example workflow

        size_t nel = 50;

        // Common parameters
        double Te = 511.0;    // keV
        double gmax = 1e3;
        double s = 2.0;

        // Test thermal distribution
        kariba::Thermal thermal(nel);
        thermal.set_temp_kev(Te);
        thermal.set_p();
        thermal.set_norm(1.0);
        thermal.set_ndens();

        CHECK(thermal.count_particles() > 0.0);
        double thermal_avg_gamma = thermal.av_gamma();
        CHECK(thermal_avg_gamma > 1.0);
        CHECK(thermal_avg_gamma < 10.0);    // Should be mildly relativistic

        // Test mixed distribution (low non-thermal fraction)
        kariba::Mixed mixed_low(nel);
        mixed_low.set_temp_kev(Te);
        mixed_low.set_p(gmax);
        mixed_low.set_plfrac(0.1);
        mixed_low.set_pspec(s);
        mixed_low.set_norm(1.0);
        mixed_low.set_ndens();

        CHECK(mixed_low.count_particles() > 0.0);

        // Test mixed distribution (high non-thermal fraction)
        kariba::Mixed mixed_high(nel);
        mixed_high.set_temp_kev(Te);
        mixed_high.set_p(gmax);
        mixed_high.set_plfrac(0.9);
        mixed_high.set_pspec(s);
        mixed_high.set_norm(1.0);
        mixed_high.set_ndens();

        CHECK(mixed_high.count_particles() > 0.0);

        // High non-thermal fraction should have higher average gamma
        double mixed_low_avg = mixed_low.av_gamma();
        double mixed_high_avg = mixed_high.av_gamma();
        CHECK(mixed_high_avg > mixed_low_avg);

        // Test kappa distribution
        kariba::Kappa kappa(nel);
        kappa.set_temp_kev(Te);
        kappa.set_p(gmax);
        kappa.set_kappa(s + 1.0);    // Kappa index convention
        kappa.set_norm(1.0);
        kappa.set_ndens();

        CHECK(kappa.count_particles() > 0.0);

        // Test broken powerlaw
        kariba::Bknpower bknpower(nel);

        // Need thermal distribution to get break momentum
        kariba::Thermal temp_for_break(nel);
        temp_for_break.set_temp_kev(Te);
        temp_for_break.set_p();
        temp_for_break.set_norm(1.0);
        temp_for_break.set_ndens();
        double pbrk = temp_for_break.av_p();

        bknpower.set_pspec1(-2.0);
        bknpower.set_pspec2(s);
        bknpower.set_p(0.1 * pbrk, pbrk, gmax);
        bknpower.set_norm(1.0);
        bknpower.set_ndens();

        CHECK(bknpower.count_particles() > 0.0);

        // Test cooling for all distributions
        double bfield = 1e3;
        double R = 1e10;
        double beta_exp = 0.1;

        CHECK_NOTHROW(mixed_low.cooling_steadystate(0.0, 1.0, bfield, R, beta_exp));
        CHECK_NOTHROW(mixed_high.cooling_steadystate(0.0, 1.0, bfield, R, beta_exp));
        CHECK_NOTHROW(kappa.cooling_steadystate(0.0, 1.0, bfield, R, beta_exp));
        CHECK_NOTHROW(bknpower.cooling_steadystate(0.0, 1.0, bfield, R, beta_exp));
    }
}

TEST_CASE("Error handling and edge cases") {
    SUBCASE("Invalid parameters") {
        kariba::Thermal thermal(10);

        // Should handle reasonable temperature range
        CHECK_NOTHROW(thermal.set_temp_kev(1.0));
        CHECK_NOTHROW(thermal.set_temp_kev(1000.0));

        // Very small or large temperatures might cause issues but shouldn't
        // crash
        CHECK_NOTHROW(thermal.set_temp_kev(0.1));
        CHECK_NOTHROW(thermal.set_temp_kev(10000.0));
    }

    SUBCASE("Array size limits") {
        // Test with minimum reasonable array size
        kariba::Thermal small(10);
        small.set_temp_kev(100.0);
        small.set_p();
        small.set_norm(1.0);
        CHECK_NOTHROW(small.set_ndens());

        // Test with larger array size
        kariba::Thermal large(500);
        large.set_temp_kev(100.0);
        large.set_p();
        large.set_norm(1.0);
        CHECK_NOTHROW(large.set_ndens());
    }

    SUBCASE("Extreme physical parameters") {
        kariba::ShSDisk disk;

        // Very small black hole
        disk.set_mbh(1.0);
        disk.set_rin(6.0 * karcst::gconst * karcst::msun / karcst::cee_cee);
        disk.set_rout(1e3 * karcst::gconst * karcst::msun / karcst::cee_cee);
        disk.set_luminosity(1e-6);
        disk.set_inclination(0.0);

        CHECK_NOTHROW(disk.disk_spectrum());

        // Very massive black hole
        disk.set_mbh(1e10);
        disk.set_rin(6.0 * karcst::gconst * 1e10 * karcst::msun / karcst::cee_cee);
        disk.set_rout(1e5 * karcst::gconst * 1e10 * karcst::msun / karcst::cee_cee);
        disk.set_luminosity(0.1);

        CHECK_NOTHROW(disk.disk_spectrum());
    }
}
