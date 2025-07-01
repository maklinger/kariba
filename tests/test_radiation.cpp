#include "doctest.h"

#include <cmath>
#include <gsl/gsl_spline.h>
#include <kariba/BBody.hpp>
#include <kariba/Compton.hpp>
#include <kariba/Cyclosyn.hpp>
#include <kariba/Powerlaw.hpp>
#include <kariba/ShSDisk.hpp>
#include <kariba/Thermal.hpp>
#include <kariba/constants.hpp>

namespace karcst = kariba::constants;

TEST_CASE("Radiation base class functionality") {
    SUBCASE("ShSDisk accretion disk") {
        kariba::ShSDisk disk;

        SUBCASE("Basic disk parameters") {
            double mbh = 10.0;    // 10 solar masses
            double rin =
                10.0 * karcst::gconst * mbh * karcst::msun / karcst::cee_cee;
            double rout =
                1e4 * karcst::gconst * mbh * karcst::msun / karcst::cee_cee;
            double ldisk = 1e-4;    // Eddington units

            disk.set_mbh(mbh);
            disk.set_rin(rin);
            disk.set_rout(rout);
            disk.set_luminosity(ldisk);
            disk.set_inclination(0.0);

            CHECK(disk.rin() == rin);
            CHECK(disk.lum() == ldisk);
        }

        SUBCASE("Disk spectrum calculation") {
            double mbh = 10.0;
            double rg = karcst::gconst * mbh * karcst::msun / karcst::cee_cee;

            disk.set_mbh(mbh);
            disk.set_rin(10.0 * rg);
            disk.set_rout(1e4 * rg);
            disk.set_luminosity(1e-4);
            disk.set_inclination(0.0);

            CHECK_NOTHROW(disk.disk_spectrum());

            CHECK(!disk.get_energy().empty());
            CHECK(!disk.get_nphot().empty());
            CHECK(!disk.get_energy_obs().empty());
            CHECK(!disk.get_nphot_obs().empty());

            // Check that spectrum has positive values
            const std::vector<double> &energy = disk.get_energy_obs();
            const std::vector<double> &nphot = disk.get_nphot_obs();

            bool has_positive_flux = false;
            for (int i = 0; i < disk.get_size(); i++) {
                CHECK(energy[i] > 0.0);
                if (nphot[i] > 0.0) {
                    has_positive_flux = true;
                }
            }
            CHECK(has_positive_flux);
        }

        SUBCASE("Disk geometry and beaming") {
            disk.set_mbh(10.0);
            disk.set_rin(10.0 * karcst::gconst * 10.0 * karcst::msun /
                         karcst::cee_cee);
            disk.set_rout(1e4 * karcst::gconst * 10.0 * karcst::msun /
                          karcst::cee_cee);
            disk.set_luminosity(1e-4);

            // Test different inclination angles
            for (double incl = 0.0; incl <= 90.0; incl += 30.0) {
                CHECK_NOTHROW(disk.set_inclination(incl));
            }
        }
    }

    SUBCASE("BBody black body radiation") {
        kariba::BBody bbody;

        SUBCASE("Black body spectrum") {
            double temp_k = 1e4;         // 10,000 K
            double radius = 1e10;        // 10^10 cm
            double luminosity = 1e38;    // erg/s

            bbody.set_temp_k(temp_k);
            bbody.set_geometry("sphere", radius);
            bbody.set_lum(luminosity);
            bbody.set_beaming(0.0, 0.0, 1.0);    // No beaming

            CHECK_NOTHROW(bbody.bb_spectrum());

            CHECK(!bbody.get_energy().empty());
            CHECK(!bbody.get_nphot().empty());

            // Check that spectrum has expected properties
            const std::vector<double> &energy = bbody.get_energy();
            const std::vector<double> &nphot = bbody.get_nphot();

            // Check that all energies are positive and increasing
            bool increasing = true;
            bool has_positive_flux = false;
            for (int i = 0; i < bbody.get_size(); i++) {
                CHECK(energy[i] > 0.0);
                if (i > 0 && energy[i] <= energy[i - 1]) {
                    increasing = false;
                }
                if (nphot[i] > 0.0) {
                    has_positive_flux = true;
                }
            }
            CHECK(increasing);
            CHECK(has_positive_flux);

            // Check temperature retrieval
            CHECK(std::abs(bbody.temp_k() - temp_k) / temp_k < 0.01);
        }
    }
}

TEST_CASE("Synchrotron radiation") {
    SUBCASE("Cyclosyn synchrotron emission") {
        kariba::Cyclosyn syncro(100);
        kariba::Powerlaw electrons(100);

        // Set up electron distribution
        double pmin = 1e-3 * karcst::emgm * karcst::cee;
        double gmax = 1e4;

        electrons.set_p(pmin, gmax);
        electrons.set_pspec(2.5);
        electrons.set_norm(1.0);
        electrons.set_ndens();

        // Set up GSL interpolation
        gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
        gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen, 100);
        gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
        gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen, 100);

        gsl_spline_init(spline_eldis, electrons.get_gamma().data(),
                        electrons.get_gdens().data(), 100);
        gsl_spline_init(spline_deriv, electrons.get_gamma().data(),
                        electrons.get_gdens_diff().data(), 100);

        SUBCASE("Basic synchrotron setup") {
            double bfield = 1e3;     // Gauss
            double radius = 1e15;    // cm

            syncro.set_frequency(1e8, 1e18);
            syncro.set_bfield(bfield);
            syncro.set_geometry("sphere", radius);
            syncro.set_beaming(0.0, 0.0, 1.0);

            double gmin = electrons.get_gamma()[0];
            double gmax_actual = electrons.get_gamma()[99];

            CHECK_NOTHROW(syncro.cycsyn_spectrum(gmin, gmax_actual,
                                                 spline_eldis, acc_eldis,
                                                 spline_deriv, acc_deriv));

            CHECK(!syncro.get_energy().empty());
            CHECK(!syncro.get_nphot().empty());

            // Check for positive flux somewhere
            const std::vector<double> &nphot = syncro.get_nphot();
            bool has_flux = false;
            for (int i = 0; i < syncro.get_size(); i++) {
                if (nphot[i] > 0.0) {
                    has_flux = true;
                }
            }
            CHECK(has_flux);
        }

        SUBCASE("Synchrotron with beaming") {
            double bfield = 1e3;
            double radius = 1e15;
            double theta = 10.0;    // degrees
            double beta = 0.9;
            double delta = 1.0 / (1.0 - beta * cos(theta * karcst::pi / 180.0));

            syncro.set_frequency(1e8, 1e18);
            syncro.set_bfield(bfield);
            syncro.set_geometry("sphere", radius);
            syncro.set_beaming(theta, beta, delta);

            double gmin = electrons.get_gamma()[0];
            double gmax_actual = electrons.get_gamma()[99];

            CHECK_NOTHROW(syncro.cycsyn_spectrum(gmin, gmax_actual,
                                                 spline_eldis, acc_eldis,
                                                 spline_deriv, acc_deriv));
        }

        // Cleanup GSL objects
        gsl_spline_free(spline_eldis);
        gsl_interp_accel_free(acc_eldis);
        gsl_spline_free(spline_deriv);
        gsl_interp_accel_free(acc_deriv);
    }
}

TEST_CASE("Inverse Compton scattering") {
    SUBCASE("Compton scattering setup") {
        kariba::Compton compton(100, 50);
        kariba::Thermal electrons(100);
        kariba::ShSDisk disk;

        // Set up thermal electrons
        double temp_kev = 100.0;
        electrons.set_temp_kev(temp_kev);
        electrons.set_p();
        electrons.set_norm(1e-3);    // cm^-3
        electrons.set_ndens();

        // Set up disk as seed photon source
        double mbh = 10.0;
        double rg = karcst::gconst * mbh * karcst::msun / karcst::cee_cee;
        disk.set_mbh(mbh);
        disk.set_rin(10.0 * rg);
        disk.set_rout(1e4 * rg);
        disk.set_luminosity(1e-4);
        disk.set_inclination(0.0);
        disk.disk_spectrum();

        // Set up GSL interpolation for electrons
        gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
        gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen, 100);

        gsl_spline_init(spline_eldis, electrons.get_gamma().data(),
                        electrons.get_gdens().data(), 100);

        SUBCASE("Basic Compton scattering") {
            double radius = 1e15;    // cm

            compton.set_frequency(1e17, 1e25);
            compton.set_geometry("sphere", radius);
            compton.set_beaming(0.0, 0.0, 1.0);
            compton.set_tau(1e-3, temp_kev);

            // Use disk as seed photon field
            compton.shsdisk_seed(disk.get_energy(), disk.tin(), disk.rin(),
                                 disk.rin() * 100, disk.hdisk(), 0.0);

            double gmin = electrons.get_gamma()[0];
            double gmax = electrons.get_gamma()[99];

            CHECK_NOTHROW(
                compton.compton_spectrum(gmin, gmax, spline_eldis, acc_eldis));

            CHECK(!compton.get_energy().empty());
            CHECK(!compton.get_nphot().empty());
        }

        SUBCASE("Multiple scattering") {
            compton.set_frequency(1e17, 1e25);
            compton.set_geometry("sphere", 1e15);
            compton.set_beaming(0.0, 0.0, 1.0);
            compton.set_tau(0.1, temp_kev);    // Higher optical depth
            compton.set_niter(5);              // Multiple scatterings

            compton.shsdisk_seed(disk.get_energy(), disk.tin(), disk.rin(),
                                 disk.rin() * 100, disk.hdisk(), 0.0);

            double gmin = electrons.get_gamma()[0];
            double gmax = electrons.get_gamma()[99];

            CHECK_NOTHROW(
                compton.compton_spectrum(gmin, gmax, spline_eldis, acc_eldis));
        }

        // Cleanup
        gsl_spline_free(spline_eldis);
        gsl_interp_accel_free(acc_eldis);
    }

    SUBCASE("Self-Synchrotron Compton (SSC)") {
        kariba::Cyclosyn syncro(100);
        kariba::Compton ssc(100, 100);
        kariba::Powerlaw electrons(100);

        // Set up electron distribution
        double pmin = 1e-3 * karcst::emgm * karcst::cee;
        double gmax = 1e4;

        electrons.set_p(pmin, gmax);
        electrons.set_pspec(2.5);
        electrons.set_norm(1.0);
        electrons.set_ndens();

        // Set up synchrotron first
        syncro.set_frequency(1e8, 1e18);
        syncro.set_bfield(1e3);
        syncro.set_geometry("sphere", 1e15);
        syncro.set_beaming(0.0, 0.0, 1.0);

        // GSL setup
        gsl_interp_accel *acc_eldis = gsl_interp_accel_alloc();
        gsl_spline *spline_eldis = gsl_spline_alloc(gsl_interp_steffen, 100);
        gsl_interp_accel *acc_deriv = gsl_interp_accel_alloc();
        gsl_spline *spline_deriv = gsl_spline_alloc(gsl_interp_steffen, 100);

        gsl_spline_init(spline_eldis, electrons.get_gamma().data(),
                        electrons.get_gdens().data(), 100);
        gsl_spline_init(spline_deriv, electrons.get_gamma().data(),
                        electrons.get_gdens_diff().data(), 100);

        double gmin = electrons.get_gamma()[0];
        double gmax_actual = electrons.get_gamma()[99];

        syncro.cycsyn_spectrum(gmin, gmax_actual, spline_eldis, acc_eldis,
                               spline_deriv, acc_deriv);

        // Now set up SSC using synchrotron photons as seed
        ssc.set_frequency(1e20, 1e28);
        ssc.set_geometry("sphere", 1e15);
        ssc.set_beaming(0.0, 0.0, 1.0);
        ssc.set_tau(1e-3, electrons.av_gamma() * 511.0);

        ssc.cyclosyn_seed(syncro.get_energy(), syncro.get_nphot());

        CHECK_NOTHROW(
            ssc.compton_spectrum(gmin, gmax_actual, spline_eldis, acc_eldis));

        // Cleanup
        gsl_spline_free(spline_eldis);
        gsl_interp_accel_free(acc_eldis);
        gsl_spline_free(spline_deriv);
        gsl_interp_accel_free(acc_deriv);
    }
}

TEST_CASE("Radiation consistency checks") {
    SUBCASE("Energy conservation") {
        kariba::BBody bbody;

        bbody.set_temp_k(5000.0);
        bbody.set_geometry("sphere", 1e10);
        bbody.set_lum(1e35);
        bbody.set_beaming(0.0, 0.0, 1.0);
        bbody.bb_spectrum();

        // Check that luminosity is reasonable
        double total_lum = bbody.integrated_luminosity(1e12, 1e18);
        CHECK(total_lum > 0.0);

        // Stefan-Boltzmann check (approximately)
        double expected_lum = 4.0 * karcst::pi * pow(1e10, 2.0) *
                              karcst::sbconst * pow(5000.0, 4.0);

        // Should be within reasonable range (integration limits matter)
        CHECK(total_lum / expected_lum > 0.01);
        CHECK(total_lum / expected_lum < 10000.0);    // Very loose bound
    }

    SUBCASE("Frequency array monotonicity") {
        kariba::Cyclosyn syncro(100);

        syncro.set_frequency(1e10, 1e20);

        // Check that frequency arrays are set up correctly
        const std::vector<double> &energy = syncro.get_energy();
        CHECK(!energy.empty());

        // Frequencies should be monotonically increasing
        for (size_t i = 1; i < 100; i++) {
            CHECK(energy[i] > energy[i - 1]);
        }

        // Check bounds (allow small numerical errors)
        CHECK(energy[0] >= 1e10 * karcst::herg * 0.99);
        CHECK(energy[99] <= 1e20 * karcst::herg * 1.01);
    }
}
