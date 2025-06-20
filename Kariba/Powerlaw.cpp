#include <fstream>
#include <iomanip>
#include <iostream>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "kariba/Particles.hpp"
#include "kariba/Powerlaw.hpp"

// Class constructor to initialize object
Powerlaw::Powerlaw(int s) {
    size = s;

    p = new double[size];
    ndens = new double[size];
    gamma = new double[size];
    gdens = new double[size];
    gdens_diff = new double[size];

    plnorm = 1.;

    mass_gr = emgm;
    mass_kev = emgm * gr_to_kev;

    for (int i = 0; i < size; i++) {
        p[i] = 0;
        ndens[i] = 0;
        gamma[i] = 0;
        gdens[i] = 0;
        gdens_diff[i] = 0;
    }
}

// Methods to set momentum/energy arrays
void Powerlaw::set_p(double min, double ucom, double bfield, double betaeff,
                     double r, double fsc) {
    pmin = min;
    pmax = max_p(ucom, bfield, betaeff, r, fsc);

    double pinc = (log10(pmax) - log10(pmin)) / (size - 1);

    for (int i = 0; i < size; i++) {
        p[i] = pow(10., log10(pmin) + i * pinc);
        gamma[i] = pow(pow(p[i] / (mass_gr * cee), 2.) + 1., 1. / 2.);
    }
}

void Powerlaw::set_p(double min, double gmax) {
    pmin = min;
    pmax = pow(pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * cee;

    double pinc = (log10(pmax) - log10(pmin)) / (size - 1);

    for (int i = 0; i < size; i++) {
        p[i] = pow(10., log10(pmin) + i * pinc);
        gamma[i] = pow(pow(p[i] / (mass_gr * cee), 2.) + 1., 1. / 2.);
    }
}

// Method to set differential electron number density from known pspec,
// normalization, and momentum array
void Powerlaw::set_ndens() {
    for (int i = 0; i < size; i++) {
        ndens[i] = plnorm * pow(p[i], -pspec) * exp(-p[i] / pmax);
    }
    initialize_gdens();
    gdens_differentiate();
}

// methods to set the slope and normalization
void Powerlaw::set_pspec(double s1) { pspec = s1; }

void Powerlaw::set_norm(double n) {
    plnorm = n * (1. - pspec) /
             (pow(p[size - 1], (1. - pspec)) - pow(p[0], (1. - pspec)));
}

// Injection function to be integrated in cooling
double injection_pl_int(double x, void *p) {
    struct injection_pl_params *params = (struct injection_pl_params *) p;
    double s = (params->s);
    double n = (params->n);
    double m = (params->m);
    double max = (params->max);

    double mom_int = pow(pow(x, 2.) - 1., 1. / 2.) * m * cee;

    return n * pow(mom_int, -s) * exp(-mom_int / max);
}

// Method to solve steady state continuity equation. NOTE: KN cross section not
// included in IC cooling
void Powerlaw::cooling_steadystate(double ucom, double n0, double bfield,
                                   double r, double betaeff) {
    double Urad = pow(bfield, 2.) / (8. * pi) + ucom;
    double pdot_ad = betaeff * cee / r;
    double pdot_rad =
        (4. * sigtom * cee * Urad) / (3. * mass_gr * pow(cee, 2.));
    double tinj = r / (cee);

    double integral, error;
    gsl_function F1;
    struct injection_pl_params params = {pspec, plnorm, mass_gr, pmax};
    F1.function = &injection_pl_int;
    F1.params = &params;

    for (int i = 0; i < size; i++) {
        if (i < size - 1) {
            gsl_integration_workspace *w1;
            w1 = gsl_integration_workspace_alloc(100);
            gsl_integration_qag(&F1, gamma[i], gamma[i + 1], 1e1, 1e1, 100, 1,
                                w1, &integral, &error);
            gsl_integration_workspace_free(w1);

            ndens[i] = (integral / tinj) /
                       (pdot_ad * p[i] / (mass_gr * cee) +
                        pdot_rad * (gamma[i] * p[i] / (mass_gr * cee)));
        } else {
            ndens[size - 1] = ndens[size - 2] *
                              pow(p[size - 1] / p[size - 2], -pspec - 1) *
                              exp(-1.);
        }
    }
    // the last bin is set by arbitrarily assuming cooled distribution; this is
    // necessary because the integral
    // above is undefined for the last bin

    // The last step requires a renormalization. The reason is that the result
    // of gsl_integration_qag strongly depends on the value of "size". Without
    // doing anything fancy, this can be fixed simply by ensuring that the total
    // integrated number of density equals n0 (which we know), and rescaling the
    // array ndens[i] by the appropriate constant.
    double renorm = count_particles() / n0;

    for (int i = 0; i < size; i++) {
        ndens[i] = ndens[i] / renorm;
    }

    initialize_gdens();
    gdens_differentiate();
}

// Method to calculate maximum momentum of non thermal particles based on
// acceleration and cooling timescales The estimate is identical to the old
// agnjet but in momentum space; see Lucchini et al. 2019 for the math of the
// old version
double Powerlaw::max_p(double ucom, double bfield, double betaeff, double r,
                       double fsc) {
    double Urad, escom, accon, syncon, b, c, gmax;
    Urad = pow(bfield, 2.) / (8. * pi) + ucom;
    escom = betaeff * cee / r;
    syncon = (4. * sigtom * Urad) / (3. * mass_gr * cee);
    accon = (3. * fsc * charg * bfield) / (4. * mass_gr * cee);

    b = escom / syncon;
    c = accon / syncon;

    gmax = (-b + pow(pow(b, 2.) + 4. * c, 1. / 2.)) / 2.;

    return pow(pow(gmax, 2.) - 1., 1. / 2.) * mass_gr * cee;
}

// Methods to set energy array for protons
void Powerlaw::set_energy(double gpmin, double fsc, double f_beta,
                          double bfield, double r_g, double z, double r,
                          int infosw, double protdens, double ntargets,
                          double Uradjet, std::string outputConfiguration,
                          std::string source) {

    double gpmax, logdgp;
    if (fsc < 1.) {
        double pp_targets = ntargets + protdens;

        ProtonTimescales(logdgp, fsc, f_beta, bfield, gpmin, gpmax, r_g, z, r,
                         infosw, pp_targets, Uradjet, outputConfiguration,
                         source);

    } else {
        isEfficient = true;
        gpmax = fsc;
        logdgp = log10(2. * gpmax / gpmin) /
                 (size - 1);    // so as to extend a bit further from γ_p,max
        check_secondary_charged_syn(bfield, gpmax);
    }

    for (int i = 0; i < size; i++) {
        gamma[i] = pow(10., (log10(gpmin) + i * logdgp));
        p[i] = sqrt(gamma[i] * gamma[i] - 1.) * mass_gr * cee;
    }
    pmin = p[0];
    pmax = sqrt(gpmax * gpmax - 1.) * mass_gr * cee;
}

void Powerlaw::ProtonTimescales(double &logdgp, double fsc, double f_beta,
                                double bfield, double gpmin, double &gpmax,
                                double r_g, double z, double r, int infosw,
                                double pp_targets, double Uradjet,
                                std::string outputConfiguration,
                                std::string source) {

    double Tacc0;    // Tacc = Ep*Tacc,0 ==> Tacc,0 = 1/(ηecB)
    double paramG,
        paramH;        // The parameters used to solve Emax equation, i.e.,
                       //  Emax^2 + paramG*Emax - paramH =0 for protons
    double betap;      // the velocity in c of the proton
    double sinel;      // the pp cross section from Kelner et al.2006
    double tescape;    // timescale of the escape especially for the low energy
                       // ones
    double EpTeV;      // energy of the proton in TeV
    double gp;         // Lorentz factor of non-thermal protons

    std::ofstream timescalesFile;

    Tacc0 = 1. / (3. * fsc / 4. * charg * cee * bfield);    // acceleration
    double Tesc = r / (f_beta * cee);                       // proton escape
    double Tpp = 1. / (sigmapp * pp_targets * cee);         // pp

    // The only photon field here is the companion's/external boosted radiation
    double Tpg0 = 1. / (9.38 * 7.4e-17 * Uradjet / (mass_gr * cee * cee));

    double Tsynp0 = 6. * pi / (sigtom * bfield * bfield) *
                    pow(mass_gr * cee, 4) /
                    (emgm * emgm * cee);    // proton synchrotron

    // The parameters used to solve Emax equation, i.e., Emax^2 + paramG*Emax -
    // paramH = 0 for protons:
    paramG = (1. / Tesc + 1. / Tpp) / (1. / Tpg0 + 1. / Tsynp0);
    paramH = (1. / Tacc0) / (1. / Tpg0 + 1. / Tsynp0);

    // Maximun proton energy in erg as calculated by losses:
    double Epmax = (-paramG + sqrt(paramG * paramG + 4. * paramH)) / 2.;

    gpmax = Epmax / (pmgm * cee * cee) + 1.;
    if (Epmax > pmgm * cee * cee) {
        isEfficient = true;
        logdgp = log10(2. * gpmax / gpmin) /
                 (size - 1);    // so as to extend a bit further from γ_p,max
    } else {
        isEfficient = false;
        logdgp =
            log10(10.) / (size - 1);    // make an array between 1 and 10 GeV
    }
    if (infosw >= 2) {
        std::string filepath = outputConfiguration +
                               "/Output/Particles/timescales_" + source +
                               ".dat";
        timescalesFile.open(filepath, std::ios::app);
        for (int i = 0; i < size; i++) {
            gp = pow(10., (log10(gpmin) + i * logdgp));
            betap = sqrt(gp * gp - 1.) / gp;
            tescape = r / (cee * betap * f_beta);

            EpTeV = pmgm * cee * cee * gp / 1.6;
            sinel = sigma_pp(EpTeV);

            // we set the correct Tpg0 in order to plot the timescales;
            Tpg0 = (gp * mass_gr * cee * cee * erg >= 4.8e14)
                       ? 1. / (9.38 * 7.4e-17 * Uradjet / (mass_gr * cee * cee))
                       : 1.e100;
            Tpp = 1. / (Kpp * mbarn * sinel * pp_targets * cee);

            timescalesFile << std::left << std::setw(13) << z / r_g
                           << std::setw(13) << gp << std::setw(13)
                           << Tacc0 * mass_gr * cee * cee * gp << std::setw(13)
                           << tescape << std::setw(13)
                           << Tsynp0 /
                                  (mass_gr * cee * cee * gp * betap * betap)
                           << std::setw(13) << Tpp << std::setw(15)
                           << Tpg0 / (mass_gr * cee * cee) / gp;
            timescalesFile << std::endl;
        }
        timescalesFile.close();
        check_secondary_charged_syn(bfield, gpmax);
    }
}

void Powerlaw::check_secondary_charged_syn(double bfield, double gpmax) {
    if ((bfield * gpmax) >
        7.8e11 *
            10.) {    // I *10 because the γ_p,max is in the exponential cutoff
        std::cout << "The pion synchrotron should be taken into account"
                  << std::endl;
        std::cout << "\tB =" << std::setw(15) << bfield << std::setw(10)
                  << " g_p,max=" << std::setw(15) << gpmax << std::endl;
    }
    if ((bfield * gpmax) > 5.6e10 * 10.) {
        std::cout << "The muon synchrotron should be taken into account"
                  << std::endl;
        std::cout << "\tB =" << std::setw(15) << bfield << std::setw(10)
                  << " g_p,max=" << std::setw(15) << gpmax << std::endl;
    }
}

double Powerlaw::sigma_pp(double Ep) {    // cross section of pp in mb (that's
                                          // why I multiply with 1.e-27 at Phie)
    double Ethres =
        1.22e-3;    // threshold energy (in TeV) for pp interactions (= 1.22GeV)
    double L = log(Ep);    // L = ln(Ep/1TeV) as definied in Kelner et al. 2006
                           // for the cross section
    double sinel;          // σ_inel in mb

    sinel = (Ep >= Ethres) ? (34.3 + 1.88 * L + 0.25 * L * L) *
                                 (1. - pow((Ethres / Ep), 4)) *
                                 (1. - pow((Ethres / Ep), 4))
                           : 1.e-80;
    return sinel;
}
double Powerlaw::set_normprot(double nprot) {
    double Epmin = mass_gr * cee * cee * gamma[0];
    double Epmax = sqrt((pmax * cee) * (pmax * cee) +
                        (mass_gr * cee * cee) * (mass_gr * cee * cee));

    // If break energy is higher than maximum energy, set normalisation for
    // uncooled distribution:
    return nprot * (1. - pspec) /
           (pow(Epmax, 1. - pspec) - pow(Epmin, 1. - pspec));
}

// Method to set differential proton number density per γ from known pspec,
// normalization and γ array
void Powerlaw::set_gdens(double r, double protdens, double nwind, double bfield,
                         double plfrac, double Uradjet) {
    /* if plfrac_p>0, namely a frac of thermal of protons accelerate*/
    // plnormprot is in #/cm3/erg/sec and gdens is in #/cm3/gamma_p
    if (isEfficient) {
        double plnormprot = set_normprot(protdens) * plfrac * cee / r;
        double Tsynp0, Tchar, Tpp;
        double betap;
        double gpmax = sqrt(pmax * pmax / (mass_gr * cee * mass_gr * cee) + 1.);
        Tsynp0 = 6. * pi / (sigtom * bfield * bfield) * pow(mass_gr * cee, 4) /
                 (emgm * emgm * cee);    // proton synchrotron
        double Tpg0 = 1. / (9.38 * 7.4e-17 * Uradjet / (mass_gr * cee * cee));
        for (int i = 0; i < size; i++) {
            betap = sqrt(gamma[i] * gamma[i] - 1.) / gamma[i];
            Tpp = 1. /
                  (nwind * sigma_pp(pmgm * cee * cee * gamma[i] / 1.6) * betap *
                   cee * mbarn);    // proton energy in TeV for the function
            Tchar = pow(cee / r + 1. / Tpp +
                            gamma[i] * mass_gr * cee * cee *
                                (1. / Tpg0 + 1. / Tsynp0),
                        -1);
            gdens[i] = plnormprot *
                       pow(gamma[i] * mass_gr * cee * cee, -pspec) *
                       exp(-gamma[i] / gpmax) * Tchar * mass_gr * cee * cee;
        }
    } else {    // in case Emax<Emin
        for (int i = 0; i < size; i++)
            gdens[i] = 1.e-100;
    }
}

// if (Lumsw == 1)
void Powerlaw::set_gdens(double r, double beta, double Ljet, double ep,
                         double pspec, double &protdens, double Urad) {
    /*Sets the normalization of the accelerated protons assuming that a fraction
    --ep--  of the jet power goes into the proton acceleration, unless Emax<Emin
    so set to zero */
    if (isEfficient) {
        double G_jet = 1. / sqrt(1. - beta * beta);    // bulk Lorentz factor
        double plnormprot;                             // in #/cm3
        double gpmax = sqrt(pmax * pmax / (mass_gr * cee * mass_gr * cee) + 1.);
        double sum = 0;
        double dx = log10(gamma[2] / gamma[1]);
        for (int i = 0; i < size; i++)
            sum += log(10.) * dx * pow(gamma[i], -pspec + 2) *
                   exp(-gamma[i] / gpmax);

        plnormprot =
            ep * Ljet /
            (mass_gr * cee * cee * pi * r * r * G_jet * beta * cee * sum);
        for (int i = 0; i < size; i++)
            gdens[i] =
                plnormprot * pow(gamma[i], -pspec) * exp(-gamma[i] / gpmax);

        sum = 0.;
        for (int i = 0; i < size; i++)
            sum += log(10.) * gdens[i] * gamma[i] * dx;
        protdens = sum;
    } else {
        for (int i = 0; i < size; i++)
            gdens[i] = 1.e-100;
        protdens = 1.e-100;
    }
}
void Powerlaw::set_gdens(double &plfrac_p, double Up, double protdens) {
    /* If mass-loading, I use the specific enthalpy to work the normalisation */
    if (isEfficient) {
        double gpmax = sqrt(pmax * pmax / (mass_gr * cee * mass_gr * cee) + 1.);
        double sum = 0;
        double dx = log10(gamma[2] / gamma[1]);
        for (int i = 0; i < size; i++)
            sum += log(10.) * pow(gamma[i], -pspec + 2.) *
                   exp(-gamma[i] / gpmax) * dx;
        double K = std::max(Up / (sum * mass_gr * cee * cee), 0.);

        sum = 0.;
        for (int i = 0; i < size; i++)
            sum += log(10.) * pow(gamma[i], -pspec + 1.) *
                   exp(-gamma[i] / gpmax) * dx;
        double n_nth = K * sum;
        plfrac_p = n_nth / protdens;

        for (int i = 0; i < size; i++)
            gdens[i] = K * pow(gamma[i], -pspec) * exp(-gamma[i] / gpmax);
    } else {
        plfrac_p = 0;
        for (int i = 0; i < size; i++)
            gdens[i] = 1.e-100;
    }
}

// simple method to check quantities.
void Powerlaw::test() {
    std::cout << "Power-law distribution;" << std::endl;
    std::cout << "pspec: " << pspec << std::endl;
    std::cout << "Array size: " << size << std::endl;
    std::cout << "Default normalization: " << plnorm << std::endl;
    std::cout << "Particle mass in grams: " << mass_gr << std::endl;
}
