#include <array>
#include <cmath>

#include <kariba/Radiation.hpp>
#include <kariba/constants.hpp>

#include "bhjet.hpp"

namespace karcst = kariba::constants;

// Velocity profile function: depending on mode either interpolate 1-d hydro jet
// velocity tables or calculate 1-d magnetically accelerated jet velocity
// profile returns x axis in either units of initial jet radius (velprof_iso and
// velprof_ad) or rg (velprof_mag) and y axis in units of gamma*beta, where
// gamma is the jet bulk lorentz factor and beta its speed in units of c For
// information on the 2 velocity profiles below see Crumley et al. 2016
void velprof_iso(gsl_spline* spline) {
    // Tabulated velocity for 1D quasi-isothermal Bernoulli eq.
    const size_t size = 54;

    std::array<double, size> gbx_vel_iso = {
        1.0,         1.00001,     1.00005,     1.00023,     1.00101,     1.00456,     1.02053,
        1.09237,     1.41567,     2.87053,     6.26251,     14.3691,     34.0825,     82.8831,
        205.572,     518.283,     1325.11,     3429.81,     8975.1,      23719.5,     63255.,
        170099.,     460962.,     1.25824e+06, 3.4578e+06,  9.5633e+06,  2.66095e+07, 7.44655e+07,
        2.09528e+08, 5.92636e+08, 1.6846e+09,  4.81148e+09, 1.38055e+10, 3.9787e+10,  1.15153e+11,
        3.34651e+11, 9.76408e+11, 2.85981e+12, 8.40728e+12, 2.4805e+13,  7.34419e+13, 2.18185e+14,
        6.50341e+14, 1.94471e+15, 5.83352e+15, 1.75523e+16, 5.29701e+16, 1.60321e+17, 4.86616e+17,
        1.48111e+18, 4.52032e+18, 1.38326e+19, 4.24394e+19, 1.0e+20};

    std::array<double, size> gby_vel_iso = {
        0.485071, 1.05031, 1.05032, 1.05039, 1.05067, 1.05193, 1.05751, 1.08105, 1.16389,
        1.35278,  1.52406, 1.68077, 1.82495, 1.95888, 2.08429, 2.20255, 2.3147,  2.42158,
        2.52386,  2.62205, 2.71662, 2.80793, 2.89628, 2.98195, 3.06516, 3.14611, 3.22497,
        3.30189,  3.37702, 3.45046, 3.52232, 3.59270, 3.66169, 3.72937, 3.79580, 3.86106,
        3.92520,  3.98827, 4.05033, 4.11143, 4.17160, 4.23089, 4.28933, 4.34696, 4.40381,
        4.45992,  4.51530, 4.56999, 4.62402, 4.67740, 4.73015, 4.78230, 4.83388, 4.87281};

    gsl_spline_init(spline, gbx_vel_iso.data(), gby_vel_iso.data(), size);
}

void velprof_ad(gsl_spline* spline) {
    // Tabluated velocity for 1D adiabatic Bernoulli eq.
    const size_t size = 54;

    std::array<double, size> gbx_vel_ad = {
        1.0,         1.01667,     1.03362,     1.08617,     1.12268,     1.16042,     1.23975,
        1.30278,     1.43863,     1.56259,     1.75429,     2.00234,     2.28546,     2.69630,
        3.34274,     3.94364,     4.73011,     5.48892,     6.16231,     7.27006,     8.43633,
        1.01187e+01, 1.19378e+01, 1.45572e+01, 1.89648e+01, 2.47070e+01, 2.96343e+01, 4.26328e+01,
        7.1172e+01,  1.1495e+02,  2.05019e+02, 4.31392e+02, 7.82228e+02, 1.30586e+03, 2.32905e+03,
        4.98237e+03, 7.05026e+03, 1.559e+04,   3.68259e+04, 1.07850e+05, 2.08929e+05, 8.10435e+05,
        1.51893e+06, 2.94251e+06, 7.93391e+06, 9.83604e+06, 1.e7,        1.e8,        1.e9,
        1.e10,       1.e11,       1.e12,       1e13,        1.e14};

    std::array<double, size> gby_vel_ad = {
        0.48507, 0.56599, 0.65850, 0.72925,  0.78367, 0.83265, 0.90340, 0.95238, 1.03401,
        1.09388, 1.17551, 1.27891, 1.360544, 1.45306, 1.55646, 1.6381,  1.70884, 1.76871,
        1.81769, 1.8666,  1.91565, 1.96462,  2.00816, 2.04626, 2.10068, 2.14966, 2.17688,
        2.22041, 2.27483, 2.30748, 2.33469,  2.35646, 2.36735, 2.37823, 2.38367, 2.38912,
        2.38912, 2.38912, 2.39456, 2.39456,  2.39456, 2.39456, 2.39456, 2.39456, 2.39456,
        2.39456, 2.39456, 2.39456, 2.39456,  2.39456, 2.39456, 2.39456, 2.39456, 2.39456};

    gsl_spline_init(spline, gbx_vel_ad.data(), gby_vel_ad.data(), size);
}

// For information on this velocity profile see Lucchini et al. 2018
void velprof_mag(jet_dynpars& dyn, gsl_spline* spline) {
    size_t size = 54;
    std::vector<double> gbx_vel_mag(size, 0.0);
    std::vector<double> gby_vel_mag(size, 0.0);
    double step;

    step = (std::log10(dyn.max) + 1. - std::log10(dyn.min)) / static_cast<double>(size - 3);
    for (size_t i = 0; i < size; i++) {
        gbx_vel_mag[i] = std::pow(10., std::log10(dyn.min) + static_cast<double>(i) * step);
        if (gbx_vel_mag[i] < dyn.h0) {
            gby_vel_mag[i] = dyn.gam0;
        } else if (gbx_vel_mag[i] < dyn.acc) {
            gby_vel_mag[i] = dyn.gam0 + ((dyn.gamf - dyn.gam0) /
                                         ((std::pow(dyn.acc, 0.5) - std::pow(dyn.h0, 0.5)))) *
                                            (std::pow(gbx_vel_mag[i], 0.5) - std::pow(dyn.h0, 0.5));
        } else {
            gby_vel_mag[i] = dyn.gamf;
        }
        gby_vel_mag[i] = std::sqrt(std::pow(gby_vel_mag[i], 2.) - 1.);
    }

    gsl_spline_init(spline, gbx_vel_mag.data(), gby_vel_mag.data(), size);
}

// Equipartition functions: calculate bfield,lepton number density,proton number
// density at the base for given jet power, jet base radius, initial speed,
// initial plasma beta, accounting for 1 or 2 jets. One of two functions is
// called depending on the base assumptions (magnetic/thermal driven jet)
void equipartition(int npsw, jet_dynpars& dyn, jet_enpars& en) {
    double eq_fac, dyn_fac;    // the two numbers that change equipartition are
                               // the jet dynamics and equipartition assumptions
    if (npsw == 0) {
        eq_fac = en.av_gamma * karcst::emerg * (1. + 1. / en.pbeta);
        dyn_fac = 2. * karcst::pi * std::pow(dyn.r0, 2.) * dyn.beta0 * dyn.gam0 * karcst::cee;
        en.lepdens = en.Nj / (eq_fac * dyn_fac);
        en.protdens = 0;
        en.bfield = std::pow(8. * karcst::pi * en.av_gamma * en.lepdens * karcst::emgm *
                                 std::pow(karcst::cee, 2.) / en.pbeta,
                             1. / 2.);
    } else if (npsw == 1) {
        eq_fac = 2. * en.av_gamma * karcst::emerg * (1. + 1. / en.pbeta);
        dyn_fac = 2. * karcst::pi * std::pow(dyn.r0, 2.) * dyn.beta0 * dyn.gam0 * karcst::cee;
        en.lepdens = en.Nj / (eq_fac * dyn_fac);
        en.protdens =
            (1. + 1. / en.pbeta) * en.av_gamma * en.lepdens * (karcst::emgm / karcst::pmgm);
        en.bfield = std::pow(8. * karcst::pi * en.av_gamma * en.lepdens * karcst::emgm *
                                 std::pow(karcst::cee, 2.) / en.pbeta,
                             1. / 2.);
    } else if (npsw == 2) {
        eq_fac = karcst::emerg * (karcst::pmgm / karcst::emgm + en.av_gamma * (1. + 1. / en.pbeta));
        dyn_fac = 2. * karcst::pi * std::pow(dyn.r0, 2.) * dyn.beta0 * dyn.gam0 * karcst::cee;
        en.lepdens = en.Nj / (eq_fac * dyn_fac);
        en.protdens = en.lepdens;
        en.bfield = std::pow(8. * karcst::pi * en.av_gamma * en.lepdens * karcst::emgm *
                                 std::pow(karcst::cee, 2.) / en.pbeta,
                             1. / 2.);
    }
    en.eta = en.lepdens / en.protdens;
    en.sig0 = std::pow(en.bfield, 2.) /
              (4. * karcst::pi * en.protdens * karcst::pmgm * std::pow(karcst::cee, 2.));
}

void equipartition(double Nj, jet_dynpars& dyn, jet_enpars& en) {
    double equip, eq_fac, dyn_fac;

    // step one: calculate proton number density from initial equipartition
    // assumptions
    en.sig0 = (1. + en.sig_acc) * dyn.gamf / dyn.gam0 - 1.;
    if (en.pbeta == 0) {
        en.eta = 1.;
        equip =
            (en.sig0 / 2.) *
            (4. / 3. + (karcst::pmgm) /
                           (en.av_gamma * karcst::emgm));    // NOTE: check energy vs lorentz factor
    } else {
        equip = 1. / en.pbeta;
        en.eta = (en.sig0 * karcst::pmgm) /
                 (en.av_gamma * karcst::emgm * (2. * equip - en.sig0 * 4. / 3.));
    }
    eq_fac = karcst::pmgm * std::pow(karcst::cee, 2.) +
             en.eta * en.av_gamma * karcst::emerg * (1. + equip);
    dyn_fac = 2. * karcst::pi * std::pow(dyn.r0, 2.) * dyn.gam0 * dyn.beta0 * karcst::cee;
    en.protdens = Nj / (eq_fac * dyn_fac);
    // step two: calculate lepton number density and magnetic field
    en.lepdens = en.eta * en.protdens;
    en.bfield = std::sqrt(8. * karcst::pi * en.lepdens * en.av_gamma * karcst::emerg * equip);
}

// Function to set up distance grid for calculations along the jet axies
void jetgrid(size_t i, grid_pars& grid, jet_dynpars& dyn, double r, double& delz, double& z) {
    double zinc, z_next;
    // note: the distance grid changes in steps of 2r up to a distance zcut, and
    // then becomes logarithmic this prevents resolution errors for the IC
    // emission near the base and allows for good resolution near up to zcut,
    // which generally includes the particle acceleration region If a pure log
    // grid was used the inverse Compton code becomes resolution-dependant
    // (Connonrs et al.2019) and the resolution near zsh becomes too low,
    // leading to bad confidence intervals for this parameter
    if (i == 0) {
        z = dyn.min;
        delz = dyn.h0;
        grid.cut = 1;
    } else if (z + 2. * r < grid.zcut) {
        z = z + delz;
        delz = 2. * r;
        grid.cut = grid.cut + 1;
    } else {
        if (i == grid.cut) {
            grid.zcut = z + delz;
        }
        zinc =
            (std::log10(dyn.max) - std::log10(grid.zcut)) / static_cast<double>(grid.nz - grid.cut);
        z = std::pow(10., std::log10(grid.zcut) + zinc * static_cast<double>(i - grid.cut));
        z_next =
            std::pow(10., std::log10(grid.zcut) + zinc * static_cast<double>(i + 1 - grid.cut));
        delz = z_next - z;
    }
}

void isojetpars(double z, jet_dynpars& dyn, jet_enpars& en, double& t, zone_pars& zone,
                gsl_spline* spline, gsl_interp_accel* acc) {

    double mj;
    double z_eval;
    double gb;
    // Note on z_eval: for the old agnjet, the velocity is given as a function
    // of jet launching point min (see Crumley et al. 2016)
    z_eval = std::log10((std::max(z - dyn.h0, 0.) + dyn.r0) / dyn.r0);

    if (z < dyn.h0) {
        gb = dyn.gam0 * dyn.beta0;
    } else {
        gb = gsl_spline_eval(spline, std::pow(10., z_eval), acc);
    }

    mj = gb / (dyn.gam0 * dyn.beta0);

    zone.gamma = std::sqrt(std::pow(gb, 2.) + 1.);
    zone.beta = std::sqrt((std::pow(zone.gamma, 2.) - 1.) / std::pow(zone.gamma, 2.));
    zone.r = dyn.r0 + std::max(z - dyn.h0, 0.) / mj;
    zone.lepdens = en.lepdens * std::pow(dyn.r0 / zone.r, 2.) / mj;

    if (z < dyn.h0) {
        zone.bfield = en.bfield * (dyn.r0 / zone.r) / std::pow(mj, 0.5);
        t = 1.;
    } else {
        zone.bfield = en.bfield * (dyn.r0 / zone.r) / std::pow(mj, 0.5 + 1. / 6.);
        t = 1. / std::pow(mj, 1. / 3.);
    }
}

void adjetpars(double z, jet_dynpars& dyn, jet_enpars& en, double& t, zone_pars& zone,
               gsl_spline* spline, gsl_interp_accel* acc) {

    double mj;
    double z_eval;
    double gb;
    // Note on z_eval: for the old agnjet, the velocity is given as a funct ion
    // of jet launching point min (see Crumley et al. 2016)
    z_eval = std::log10((std::max(z - dyn.h0, 0.) + dyn.r0) / dyn.r0);

    if (z < dyn.h0) {
        gb = dyn.gam0 * dyn.beta0;
    } else {
        gb = gsl_spline_eval(spline, std::pow(10., z_eval), acc);
    }

    mj = gb / (dyn.gam0 * dyn.beta0);

    zone.gamma = std::sqrt(std::pow(gb, 2.) + 1.);
    zone.beta = std::sqrt((std::pow(zone.gamma, 2.) - 1.) / std::pow(zone.gamma, 2.));

    zone.r = dyn.r0 + std::max(z - dyn.h0, 0.) / mj;
    zone.lepdens = en.lepdens * std::pow(dyn.r0 / zone.r, 2.) / mj;

    if (z < dyn.h0) {
        zone.bfield = en.bfield * (dyn.r0 / zone.r) / std::pow(mj, 0.5);
        t = 1.;
    } else {
        zone.bfield = en.bfield * (dyn.r0 / zone.r) / std::pow(mj, 0.5 + 1. / 6.);
        t = 1. / std::pow(mj, 1. / 3.) / std::pow(zone.r / dyn.r0, 2. / 3.);
    }
}

void bljetpars(double z, jet_dynpars& dyn, jet_enpars& en, double& t, zone_pars& zone,
               gsl_spline* spline, gsl_interp_accel* acc) {

    double mj, theta, theta_acc, n_acc, b_acc, gb, r_acc;
    double gb0 = dyn.gam0 * dyn.beta0;
    double gbf = std::sqrt(std::pow(dyn.gamf, 2.) - 1.);

    if (z < dyn.h0) {
        gb = gb0;
    } else if (z < dyn.acc) {
        gb = gsl_spline_eval(spline, z, acc);
    } else {
        gb = gbf;
    }

    mj = gb / gb0;

    zone.gamma = std::sqrt(std::pow(gb, 2.) + 1.);
    zone.beta = std::sqrt((std::pow(zone.gamma, 2.) - 1.) / std::pow(zone.gamma, 2.));

    theta = 0.15 / zone.gamma;
    zone.r = dyn.r0 + std::max(z - dyn.h0, 0.) * tan(theta);
    theta_acc = 0.15 / dyn.gamf;
    r_acc = dyn.r0 + (dyn.acc - dyn.h0) * tan(theta_acc);
    n_acc = en.lepdens * std::pow(dyn.r0 / r_acc, 2.) * (gb0 / gbf);
    zone.lepdens = en.lepdens * std::pow(dyn.r0 / zone.r, 2.) / mj;

    if (z < dyn.h0) {
        b_profile(zone.gamma, zone.lepdens, dyn, en, zone.bfield);
        t = 1.;
    } else if (z < dyn.acc) {
        t = 1.;
        b_profile(zone.gamma, zone.lepdens, dyn, en, zone.bfield);
    } else {
        t = 1.;
        b_profile(zone.gamma, n_acc, dyn, en, b_acc);
        zone.bfield = b_acc * (dyn.acc / z);
    }
}

void b_profile(double gam, double n, jet_dynpars& dyn, jet_enpars& en, double& field) {
    double sigma, w;

    // Only accounts for the injected distribution, not for extra accelerated
    // particles (Lucchini et al. 2018) This should not introduce any errors as
    // long as the average Lorentz factor of the electrons is below ~a few 10^2
    // and/or the pair content of the jet is limited
    w = 4. / 3. * en.av_gamma * n * karcst::emgm * std::pow(karcst::cee, 2.);
    sigma = (dyn.gam0 / gam) * (1. + en.sig0) - 1.;
    field = std::sqrt(sigma * 4. * karcst::pi *
                      (n / en.eta * karcst::pmgm * std::pow(karcst::cee, 2.) + w));
}

// This function is used to set up the external AGN photon fields (torus, BLR)
// as a function of accretion rate See Ghisellini and Tavecchio 2009 for details
void agn_photons_init(double lum, double f1, double f2, com_pars& agn_com) {
    agn_com.rblr = 1.e17 * std::pow(lum / 1.e45, 1. / 2.);
    agn_com.ublr =
        (17. / 12.) * f1 * lum / (4. * karcst::pi * std::pow(agn_com.rblr, 2.) * karcst::cee);
    agn_com.tblr = 10.2e-3;    // temperature of the BLR photons in kev
    agn_com.lblr =
        (12. / 17.) * 4. * karcst::pi * std::pow(agn_com.rblr, 2.) * karcst::cee * agn_com.ublr;

    agn_com.rdt = 2.5e18 * std::pow(lum / 1.e45, 1. / 2.);
    agn_com.udt = f2 * lum / (4. * karcst::pi * std::pow(agn_com.rdt, 2.) * karcst::cee);
    agn_com.tdt = 370.;    // temperature of the torus photons in Kelvin
    agn_com.ldt = 4. * karcst::pi * std::pow(agn_com.rdt, 2.) * karcst::cee * agn_com.udt;
}

// This function calculates the energy densities for the standard AGN photon
// fields (BLR and torus), taking beaming/debeaming into account. It sets the
// the total energy density in agn_com (used for cooling and for the IC
// estimate), as well as the individual energy densities, used in the actual IC
// calculations. The reason the last two are not part of agn_com is that
// agn_com.ublr and agn_com.udt are calculated in the observer frame while
// ublr_zone and udt_zone are in the jet comoving frame. The reason why ublr/udt
// are multiplied by delta^2 instead of gamma^2 as in Dermer, 1995. Because of
// this, when we compute urad_total (which goes to the cooling term) we
// sometimes need to convert again with a factor gamma/delta^2
void zone_agn_phfields(double z, zone_pars& zone, double& ublr_zone, double& udt_zone,
                       com_pars& agn_com) {
    // these are to calculate the conversion factor to account for z>zblr/zdt,
    // when photons are deboosted. See Ghisellini&Tavecchio 2009
    double blr_conv, dt_conv, mu_blr1, mu_blr2, mu_dt1, mu_dt2, fr;
    fr = 12. / 17.;    // the only reason this variable is here is to save space
                       // in the code

    if (z < agn_com.rblr) {    // for low z both photon fields are boosted
        ublr_zone = std::pow(zone.delta, 2.) * agn_com.ublr;
        udt_zone = std::pow(zone.delta, 2.) * agn_com.udt;
        agn_com.urad_total = std::pow(zone.gamma, 2.) * (agn_com.udt + agn_com.ublr);
    } else if (z < 3. * agn_com.rblr) {    // BLR starts being slightly deboosted
        double blr_boost, blr_interp;
        mu_blr1 = std::pow(10. / 9., -1. / 2.);
        mu_blr2 = std::pow(8. / 9., 1. / 2.);
        blr_conv = 2. * std::pow(1. - zone.beta * mu_blr1, 3.) -
                   std::pow(1. - zone.beta * mu_blr2, 3.) - std::pow(1. - zone.beta, 3.);
        blr_boost = std::pow(zone.delta, 2.) * agn_com.ublr;
        blr_interp = std::pow(zone.delta, 2.) * agn_com.ublr * blr_conv / (3. * zone.beta);

        ublr_zone = blr_boost + (z - agn_com.rblr) / (2. * agn_com.rblr) * (blr_interp - blr_boost);
        udt_zone = std::pow(zone.delta, 2.) * agn_com.udt;
        agn_com.urad_total = std::pow(zone.gamma, 2.) *
                             (agn_com.udt + ublr_zone * std::pow(zone.gamma / zone.delta, 2.));
    } else if (z < agn_com.rdt) {    // for intermediate z BLR is deboosted,
                                     // torus is boosted
        mu_blr1 = std::pow(1. + std::pow(agn_com.rblr / z, 2.), -1. / 2.);    //
        mu_blr2 = std::pow(1. - std::pow(agn_com.rblr / z, 2.), 1. / 2.);
        blr_conv = 2. * std::pow(1. - zone.beta * mu_blr1, 3.) -
                   std::pow(1. - zone.beta * mu_blr2, 3.) - std::pow(1. - zone.beta, 3.);

        ublr_zone = std::pow(zone.delta, 2.) * agn_com.ublr * blr_conv / (3. * zone.beta);
        udt_zone = std::pow(zone.delta, 2.) * agn_com.udt;
        agn_com.urad_total = std::pow(zone.gamma, 2.) *
                             (agn_com.udt + fr * agn_com.ublr * blr_conv / (3. * zone.beta));
    } else if (z < 3. * agn_com.rdt) {    // BLR fully deboosted, DT slightly deboosted
        mu_blr1 = std::pow(1. + std::pow(agn_com.rblr / z, 2.), -1. / 2.);
        mu_blr2 = std::pow(1. - std::pow(agn_com.rblr / z, 2.), 1. / 2.);
        blr_conv = 2. * std::pow(1. - zone.beta * mu_blr1, 3.) -
                   std::pow(1. - zone.beta * mu_blr2, 3.) - std::pow(1. - zone.beta, 3.);

        ublr_zone = std::pow(zone.delta, 2.) * agn_com.ublr * blr_conv / (3. * zone.beta);

        double dt_boost, dt_interp;
        mu_dt1 = std::pow(10. / 9., -1. / 2.);
        mu_dt2 = std::pow(8. / 9., 1. / 2.);
        dt_conv = 2. * std::pow(1. - zone.beta * mu_dt1, 3.) -
                  std::pow(1. - zone.beta * mu_dt2, 3.) - std::pow(1. - zone.beta, 3.);
        dt_boost = std::pow(zone.delta, 2.) * agn_com.udt;
        dt_interp = std::pow(zone.delta, 2.) * agn_com.udt * dt_conv / (3. * zone.beta);

        udt_zone = dt_boost + (z - agn_com.rdt) / (2. * agn_com.rdt) * (dt_interp - dt_boost);
        agn_com.urad_total =
            std::pow(zone.gamma, 2.) * (udt_zone * std::pow(zone.gamma / zone.delta, 2.) +
                                        fr * agn_com.ublr * blr_conv / (3. * zone.beta));

    } else {    // for high z both photon fields are deboosted
        mu_blr1 = std::pow(1. + std::pow(agn_com.rblr / z, 2.), -1. / 2.);
        mu_blr2 = std::pow(1. - std::pow(agn_com.rblr / z, 2.), 1. / 2.);
        blr_conv = 2. * std::pow(1. - zone.beta * mu_blr1, 3.) -
                   std::pow(1. - zone.beta * mu_blr2, 3.) - std::pow(1. - zone.beta, 3.);

        mu_dt1 = std::pow(1. + std::pow(agn_com.rdt / z, 2.), -1. / 2.);
        mu_dt2 = std::pow(1. - std::pow(agn_com.rdt / z, 2.), 1. / 2.);
        dt_conv = 2. * std::pow(1. - zone.beta * mu_dt1, 3.) -
                  std::pow(1. - zone.beta * mu_dt2, 3.) - std::pow(1. - zone.beta, 3.);

        ublr_zone = std::pow(zone.delta, 2.) * agn_com.ublr * blr_conv / (3. * zone.beta);
        udt_zone = std::pow(zone.delta, 2.) * agn_com.udt * dt_conv / (3. * zone.beta);
        agn_com.urad_total = std::pow(zone.gamma, 2.) *
                             (agn_com.udt * dt_conv + fr * agn_com.ublr * blr_conv) /
                             (3. * zone.beta);
    }
}
