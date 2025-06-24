#pragma once

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>


namespace kariba {

namespace constants {
    const double kpc = 1e3 * GSL_CONST_CGSM_PARSEC;
    const double cee = GSL_CONST_CGSM_SPEED_OF_LIGHT;
    const double cee_cee = cee * cee;
    const double emgm = GSL_CONST_CGSM_MASS_ELECTRON;
    const double pmgm = GSL_CONST_CGSM_MASS_PROTON;
    const double kboltz = GSL_CONST_CGSM_BOLTZMANN;
    const double kboltz_kev2erg = 1.6022e-9;    // Boltzman constant in keV/erg
    const double gr_to_kev = 5.6095883571872e+29;
    const double me_kev = 511.0;
    const double emerg =
        GSL_CONST_CGSM_MASS_ELECTRON * pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.0);
    const double pi = M_PI;
    const double charg = 4.8e-10;
    const double sigtom = GSL_CONST_CGSM_THOMSON_CROSS_SECTION;
    const double herg = GSL_CONST_CGSM_PLANCKS_CONSTANT_H;
    const double hkev = GSL_CONST_CGSM_PLANCKS_CONSTANT_H * 6.2415e8;
    const double mjy = 1.e-26;
    const double re0 = 2.81794e-13;
    const double gconst = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT;
    const double sbconst = GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT;
    const double aconst = 7.56e-15;
    const double msun = GSL_CONST_CGSM_SOLAR_MASS;
    const double erg = 6.24e11;               // 1 erg = 6.24e11 eV
    const double mprotTeV = 938.272046e-6;    // mass of proton in TeV/c^2
    const double mpionTeV = 139.57e-6;        // mass of pion in TeV/c^2
    const double Kpp = 0.5;    // Inelasticity Kpp. Here is considered constant.
    const double Kpi =
        0.17;    // fraction of E_kinetic of proton transferred to neutrinos
    const double hbar = herg / (2.0 * pi);    // h bar
    const double barn = 1.0e-24;
    const double mbarn = 1.e-3 * barn;
    const double sigmapp = 3.43e-26;    // pp cross section in cm2

}    // namespace constants

}    // namespace constants

