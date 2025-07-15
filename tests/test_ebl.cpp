#include "doctest.h"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include <kariba/EBL.hpp>

using kariba::ebl_atten_gil;

const double EPS = 1.0e-5;

TEST_CASE("EBL correction") {
    size_t ne = 50;
    std::vector<double> energy(ne, 0.0);
    std::vector<double> luminosity(ne, 1e40);
    double total_lum = 0.0;
    // Set of energies that cover 1e-3 to 100 TeV, but slightly
    // offset from integer values
    for (size_t i = 0; i < ne; i++) {
        energy[i] = -2.9 + 4.8 * static_cast<double>(i) / 50.0;
        energy[i] = std::pow(10, energy[i]);
    }

    SUBCASE("Base verification") {
        // Verify that the total sum is simply 50 * 1e40 = 5e41.
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        CHECK(total_lum == doctest::Approx(5e41).epsilon(EPS));
    }

    SUBCASE("Lower redshift limit") {
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        ebl_atten_gil(energy, luminosity, 0);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        // Luminosities have not changed
        CHECK(total_lum == doctest::Approx(5e41).epsilon(EPS));
    }

    SUBCASE("Upper redshift limit") {
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        ebl_atten_gil(energy, luminosity, 9.01);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        // All luminosities are 0, except for the lowest energy
        CHECK(total_lum == doctest::Approx(1e40).epsilon(EPS));
    }

    SUBCASE("Basic functionality") {
        // Beyond redshift 3, the model or interpolation is
        // problematic, if not simply faulty. We skip those for now
        // Instead, we test a small range of redshifts between 0.01
        // and 3
        const size_t nz = 9;
        std::array<double, nz> total_lums = {4.38479e+41, 4.09303e+41, 3.72399e+41,
                                             3.3029e+41,  2.90535e+41, 2.5852e+41,
                                             2.28577e+41, 1.94457e+41, 9.87031e+40};
        double z = 0.01;
        for (size_t j = 0; j < nz; z *= 2, ++j) {
            // Reset luminosities
            std::fill(luminosity.begin(), luminosity.end(), 1e40);
            ebl_atten_gil(energy, luminosity, z);
            total_lum = std::reduce(luminosity.begin(), luminosity.end());
            CHECK(total_lum == doctest::Approx(total_lums[j]).epsilon(EPS));
        }
    }

    SUBCASE("Minimal luminosity") {
        // Verify no changes for very low luminosities, below the
        // minimum luminosity of 10
        const double min_lum = 10.0;          // static const inside EBL
        const double lum = 0.99 * min_lum;    // below minimum lum.
        const double total = static_cast<double>(ne) * lum;
        std::fill(luminosity.begin(), luminosity.end(), lum);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        CHECK(total_lum == doctest::Approx(total).epsilon(EPS));
        ebl_atten_gil(energy, luminosity, 1.0);    // any redshift can work
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        // Luminosities have not changed
        CHECK(total_lum == doctest::Approx(total).epsilon(EPS));
    }

    SUBCASE("Infinity in input luminosity") {
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        luminosity[20] = std::numeric_limits<double>::infinity();
        ebl_atten_gil(energy, luminosity, 1.0);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        CHECK(std::isinf(total_lum));
    }

    SUBCASE("NaN in input luminosity") {
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        luminosity[20] = std::nan("0");
        ebl_atten_gil(energy, luminosity, 1.0);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        CHECK(std::isnan(total_lum));
    }

    SUBCASE("Negative input energy") {
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        for (auto &e : energy) {
            e = -e;
        }
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        ebl_atten_gil(energy, luminosity, 1.0);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        // No change, since negative energies are simply outside of
        // the valid range
        CHECK(total_lum == doctest::Approx(5e41).epsilon(EPS));
    }

    SUBCASE("Far too high energy") {
        for (auto &e : energy) {
            e = -1e9 * e;
        }
        std::fill(luminosity.begin(), luminosity.end(), 1e40);
        ebl_atten_gil(energy, luminosity, 1.0);
        total_lum = std::reduce(luminosity.begin(), luminosity.end());
        // No change, since the energies are simply outside of
        // the valid range
        CHECK(total_lum == doctest::Approx(5e41).epsilon(EPS));
    }

    // Verify the routine still works out if the input data is really short
    SUBCASE("2-element input") {
        ne = 2;
        energy.resize(ne);
        luminosity.resize(ne);
        energy[0] = 0.005;
        energy[1] = 0.50;
        luminosity[0] = 1e40;
        luminosity[1] = 1e40;
        ebl_atten_gil(energy, luminosity, 1.0);
        // No attenuation at low energies
        CHECK(luminosity[0] == doctest::Approx(1e40).epsilon(EPS));
        // Significant attenuation at high energies
        CHECK(luminosity[1] == doctest::Approx(1.50292e+37).epsilon(EPS));
    }

    SUBCASE("1-element input") {
        ne = 1;
        energy.resize(ne);
        luminosity.resize(ne);
        energy[0] = 1;
        luminosity[0] = 1e40;
        ebl_atten_gil(energy, luminosity, 1.0);

        CHECK(luminosity[0] == doctest::Approx(6.85603e+34).epsilon(EPS));
    }

    SUBCASE("empty (0-size) input") {
        ne = 1;
        energy.resize(0);
        luminosity.resize(0);
        ebl_atten_gil(energy, luminosity, 1.0);
    }
}
