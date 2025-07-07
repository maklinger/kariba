// #define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <kariba/Bknpower.hpp>
#include <kariba/Thermal.hpp>
#include <kariba/constants.hpp>

namespace karcst = kariba::constants;

TEST_CASE("Bknpower tests") {
    SUBCASE("Simple Bknpower usage") {

        size_t nel = 50;
        double gmax = 1e3;
        double s = 2.0;
        double pbreak = 8.65822e-17;
        double bfield = 1e3;
        double R = 1e10;
        double beta_exp = 0.1;

        // Test broken powerlaw
        kariba::Bknpower bknpower(nel);

        bknpower.set_pspec1(-2.0);
        bknpower.set_pspec2(s);
        bknpower.set_p(0.1 * pbreak, pbreak, gmax);
        bknpower.set_norm(1.0);
        bknpower.set_ndens();

        CHECK(bknpower.count_particles() > 0.0);

        // Test cooling
        CHECK_NOTHROW(bknpower.cooling_steadystate(0.0, 1.0, bfield, R, beta_exp));
    }
}
