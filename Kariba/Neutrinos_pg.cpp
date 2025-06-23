#include <iostream>

#include "kariba/Neutrinos_pg.hpp"
#include "kariba/Radiation.hpp"

namespace kariba {

#ifndef PHOTOMESON_TABLES
#define PHOTOMESON_TABLES
//----------------------- electrons -----------------------//
static const double etaeTable[10] = {3.0, 4.0, 5.0,  6.0,  7.0,
                                     8.0, 9.0, 10.0, 30.0, 100.0};
static const double seTable[10] = {0.658, 0.348, 0.286, 0.256, 0.258,
                                   0.220, 0.217, 0.192, 0.125, 0.0507};
static const double deltaeTable[10] = {3.09, 2.81, 2.39, 2.27, 2.13,
                                       2.20, 2.13, 2.19, 2.27, 2.63};
static const double BetaeTable[10] = {6.43e-19, 9.91e-18, 1.24e-16, 2.67e-16,
                                      3.50e-16, 4.03e-16, 4.48e-16, 4.78e-16,
                                      1.64e-15, 4.52e-15};

//----------------------- positrons -----------------------//
static const double etaposTable[20] = {1.1, 1.2, 1.3, 1.4,  1.5, 1.6,  1.7,
                                       1.8, 1.9, 2.0, 3.0,  4.0, 5.0,  6.0,
                                       7.0, 8.0, 9.0, 10.0, 30., 100.0};
static const double sposTable[20] = {0.367,  0.282,  0.260,  0.239,  0.224,
                                     0.207,  0.198,  0.193,  0.187,  0.181,
                                     0.122,  0.106,  0.0983, 0.0875, 0.0830,
                                     0.0783, 0.0735, 0.0644, 0.0333, 0.0224};
static const double deltaposTable[20] = {
    3.12, 2.96, 2.83, 2.76, 2.69, 2.66, 2.56, 2.52, 2.49, 2.48,
    2.50, 2.46, 2.46, 2.44, 2.44, 2.45, 2.5,  2.77, 2.86};
static const double BetaposTable[20] = {
    8.09e-19, 7.70e-18, 2.05e-17, 3.66e-17, 5.48e-17, 7.39e-17, 9.52e-17,
    1.20e-16, 1.47e-16, 1.75e-16, 3.31e-16, 4.16e-16, 5.57e-16, 6.78e-16,
    7.65e-16, 8.52e-16, 9.17e-16, 9.57e-16, 3.07e-15, 1.58e-14};

//----------------------- muon neutrinos -----------------------//
static const double eta_muonTable[20] = {1.1, 1.2, 1.3, 1.4,  1.5, 1.6,  1.7,
                                         1.8, 1.9, 2.0, 3.0,  4.0, 5.0,  6.0,
                                         7.0, 8.0, 9.0, 10.0, 30., 100.0};
static const double s_muonTable[20] = {
    0.0,   0.0778, 0.242, 0.377, 0.440, 0.450, 0.461, 0.451, 0.464,  0.446,
    0.366, 0.249,  0.204, 0.174, 0.156, 0.140, 0.121, 0.107, 0.0705, 0.0463};
static const double delta_muonTable[20] = {
    0.0,  0.306, 0.792, 1.09, 1.06, 0.953, 0.956, 0.922, 0.912, 0.940,
    1.49, 2.03,  2.18,  2.24, 2.28, 2.32,  2.39,  2.46,  2.53,  2.62};
static const double Beta_muonTable[20] = {
    1.08e-18, 9.91e-18, 2.47e-17, 4.43e-17, 6.70e-17, 9.04e-17, 1.18e-16,
    1.32e-16, 1.77e-16, 2.11e-16, 3.83e-16, 5.09e-16, 7.26e-16, 9.26e-16,
    1.07e-15, 1.19e-15, 1.29e-15, 1.40e-15, 5.65e-15, 3.01e-14};

//----------------------- muon antineutrinos -----------------------//
static const double eta_antimuonTable[20] = {
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,  1.9, 2.0,
    3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30., 100.0};
static const double s_antimuonTable[20] = {
    0.365,  0.287,  0.250,  0.238,  0.220,  0.206,  0.197,
    0.193,  0.187,  0.178,  0.123,  0.106,  0.0944, 0.0829,
    0.0801, 0.0752, 0.0680, 0.0615, 0.0361, 0.0228};
static const double delta_antimuonTable[20] = {
    3.09, 2.96, 2.89, 2.76, 2.71, 2.67, 2.62, 2.56, 2.52, 2.51,
    2.48, 2.56, 2.57, 2.58, 2.54, 2.53, 2.56, 2.60, 2.78, 2.88};
static const double Beta_atnimuonTable[20] = {
    8.09e-19, 7.70e-18, 1.99e-17, 3.62e-17, 5.39e-17, 7.39e-17, 9.48e-17,
    1.20e-16, 1.47e-16, 1.74e-16, 3.38e-16, 5.17e-16, 7.61e-16, 9.57e-16,
    1.11e-15, 1.25e-15, 1.36e-15, 1.46e-15, 5.87e-15, 3.10e-14};

//----------------------- electron neutrinos -----------------------//
static const double eta_electronTable[20] = {
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,  1.9, 2.0,
    3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30., 100.0};
static const double s_electronTable[20] = {
    0.768, 0.569, 0.491, 0.395,  0.310,  0.323,  0.305, 0.285,  0.270,  0.259,
    0.158, 0.129, 0.113, 0.0996, 0.0921, 0.0861, .0800, 0.0723, 0.0411, 0.0283};
static const double delta_electronTable[20] = {
    2.49, 2.35, 2.41, 2.45, 2.45, 2.43, 2.40, 2.39, 2.37, 2.35,
    2.42, 2.46, 2.45, 2.46, 2.46, 2.45, 2.47, 2.51, 2.70, 2.77};
static const double Beta_electronTable[20] = {
    9.43e-19, 9.22e-18, 2.35e-17, 4.20e-17, 6.26e-17, 8.57e-17, 1.13e-16,
    1.39e-16, 1.70e-16, 2.05e-16, 3.81e-16, 4.74e-16, 6.30e-16, 7.65e-16,
    8.61e-16, 9.61e-16, 1.03e-15, 1.10e-15, 3.55e-15, 1.84e-14};

//----------------------- electron antineutrinos -----------------------//
static const double eta_antielectronTable[10] = {3.0, 4.0, 5.0,  6.0,  7.0,
                                                 8.0, 9.0, 10.0, 30.0, 100.0};
static const double s_antielectronTable[10] = {
    0.985, 0.378, 0.310, 0.327, 0.308, 0.292, 0.260, 0.233, 0.135, 0.077};
static const double delta_antielectronTable[10] = {
    2.63, 2.98, 2.31, 2.11, 2.03, 1.98, 2.02, 2.07, 2.24, 2.40};
static const double Beta_antielectronTable[10] = {
    6.61e-19, 9.74e-18, 1.34e-16, 2.91e-16, 3.81e-16,
    4.48e-16, 4.83e-16, 5.13e-16, 1.75e-15, 5.48e-15};

#endif

Neutrinos_pg::~Neutrinos_pg() {
    delete[] en_phot;
    delete[] num_phot;
    delete[] en_phot_obs;
    delete[] num_phot_obs;
}

Neutrinos_pg::Neutrinos_pg(int s1, double Emin, double Emax) {

    size = s1;

    en_phot = new double[size];
    num_phot = new double[size];
    en_phot_obs = new double[2 * size];
    num_phot_obs = new double[2 * size];

    double Einc = log10(Emax / Emin) / (size - 1);
    for (int i = 0; i < size; i++) {
        en_phot[i] = pow(10., log10(Emin) + i * Einc);
        en_phot_obs[i] = en_phot[i];
        en_phot_obs[i + size] = en_phot[i];
        num_phot[i] = 0;
        num_phot_obs[i] = 0;
        num_phot_obs[i + size] = 0;
    }
}
//************************************************************************************************************
void Neutrinos_pg::set_neutrinos(double gp_min, double gp_max,
                                 gsl_interp_accel *acc_Jp,
                                 gsl_spline *spline_Jp, double *en_perseg,
                                 double *lum_perseg, int nphot,
                                 std::string outputConfiguration,
                                 std::string flavor, int infosw,
                                 std::string source) {

    std::ofstream PhotopionFile;    // for plotting
    if (infosw >= 2) {
        std::string filepath;
        if (source.compare("JET") == 0) {
            filepath =
                outputConfiguration + "/Output/Neutrinos/" + flavor + "_pg.dat";
        } else {
            std::cerr << "Wrong source; cannot be " << source
                      << " but rather JET!" << std::endl;
            exit(1);
        }
        if (not(flavor.compare("electrons") == 0 ||
                flavor.compare("positrons") == 0)) {
            PhotopionFile.open(filepath, std::ios::app);
        }
    }
    int N = 10;
    double Epion = 139.6e6 / constants::erg;    // rest mass of pion in erg
    double eta, deta;    // eta parameter: η = 4εE_p/(m_p^2*c^4) and its step
    double eta_zero = 0.313;    // eq 16 from Kelner & Aharonian 08
    double Hfunction;           // the quantity that we intergrate for every eta
    double dNdEv;               // spectrum of products in #/erg/cm3/sec
    double eta_max = 99.99;     // max η
    double eta_min = 1.10;      // min η
    double nu_min =
        en_perseg[0] / constants::herg;    // the min freq of photon targets
    double nu_max = en_perseg[nphot - 1] /
                    constants::herg;    // the max freq of photon targets
    double *freq =
        new double[nphot];    // frequency of photons per segment in Hz
    double *Uphot =
        new double[nphot];    // diff energy density per segment in #/cm3/erg

    if (flavor.compare("electrons") == 0 || flavor.compare("antielectron") == 0)
        eta_min = 3.001;
    if (flavor.compare("gamma_rays") == 0)
        Epion = 137.5e6 / constants::erg;

    for (int k = 0; k < nphot; k++) {
        freq[k] = en_perseg[k] / constants::herg;    // Hz from erg
        Uphot[k] = lum_perseg[k] * (r / constants::cee /
                                    (constants::herg * constants::herg *
                                     freq[k] * vol));    // #/cm3/erg
    }

    // Interpolation for jet photon distribution
    gsl_interp_accel *acc_ng = gsl_interp_accel_alloc();
    gsl_spline *spline_ng = gsl_spline_alloc(gsl_interp_steffen, nphot);
    gsl_spline_init(spline_ng, freq, Uphot, nphot);

    deta = log10(eta_max / eta_min) / (N - 1);
#pragma omp parallel for private(                                              \
        eta, Hfunction, dNdEv)    // possibly lost: 9,424 bytes in 31 blocks
    for (int i = 0; i < size;
         i++) {    // for every neutrino of energy en_phot[i](erg)
        double sum;
        double Ev = en_phot[i];
        if (Ev > Epion && Ev <= gp_max * constants::pmgm * constants::cee *
                                    constants::cee) {    // in erg
            sum = 0.0;
            gsl_integration_workspace *w1 =
                gsl_integration_workspace_alloc(100);
            double result1, error1;
            gsl_function F1;
            for (int j = 0; j < N; j++) {    // eq 69 from KA08
                eta = eta_zero * (pow(10., log10(eta_min) + j * deta));
                Heta_params F1params = {eta,    eta_zero,  Ev,     gp_min,
                                        gp_max, spline_Jp, acc_Jp, flavor,
                                        acc_ng, spline_ng, nu_min, nu_max};
                F1.function = &Heta;
                F1.params = &F1params;
                double max = log10(Ev / (gp_min * constants::pmgm *
                                         constants::cee * constants::cee));
                double min = log10(Ev / (gp_max * constants::pmgm *
                                         constants::cee * constants::cee));
                gsl_integration_qag(&F1, min, max, 1e0, 1e0, 100, 1, w1,
                                    &result1, &error1);
                // Have to increase to 3 to get a good shape without arificial
                // features
                Hfunction =
                    pow(constants::pmgm * constants::cee * constants::cee, 2) /
                    4. * result1;    // eq 70 from KA08
                sum += Hfunction * deta * eta * log(10.);
            }
            dNdEv = sum;    // in #/erg/cm3/sec
            gsl_integration_workspace_free(w1);
        } else {
            dNdEv = 1.e-100;    // in #/erg/cm3/sec
        }
        num_phot[i] =
            dNdEv * constants::herg * en_phot[i] * vol;    // erg/sec/Hz
        en_phot_obs[i] = en_phot[i] * dopfac;              //*dopfac;
        num_phot_obs[i] =
            num_phot[i] * pow(dopfac, dopnum);    //*dopfac*dopfac;
                                                  ////L'_v' -> L_v
    }

    if ((infosw >= 2) && !(flavor.compare("electrons") == 0 ||
                           flavor.compare("positrons") == 0)) {
        for (int i = 0; i < size; i++) {
            PhotopionFile << std::left << std::setw(15) << en_phot[i]
                          << std::setw(25)
                          << num_phot[i] / (constants::herg * en_phot[i] * vol)
                          << std::setw(25)
                          << num_phot[i] / (constants::herg * en_phot[i])
                          << std::endl;
        }
    }

    if (infosw >= 2)
        PhotopionFile.close();
    gsl_spline_free(spline_ng);
    gsl_interp_accel_free(acc_ng);
    delete[] freq;
    delete[] Uphot;
}
double Heta(double x, void *p) {
    // eq 70 from KA08 for pairs and writen as {0< x=Ee/Ep <1} Ee/Epmax <
    // x=Ee/Ep <Ee/Epmin
    Heta_params *params = (Heta_params *) p;
    double eta = (params->eta);
    double eta_zero = (params->eta_zero);
    double E = (params->E);
    double gp_min = (params->gp_min);
    double gp_max = (params->gp_max);
    gsl_spline *spline_Jp = (params->spline_Jp);
    gsl_interp_accel *acc_Jp = (params->acc_Jp);
    std::string product = (params->product);
    gsl_interp_accel *acc_ng = (params->acc_ng);
    gsl_spline *spline_ng = (params->spline_ng);
    double nu_min = (params->nu_min);
    double nu_max = (params->nu_max);

    double fp;         // number density of accelerated protons in #/cm3/erg
    double fph;        // number density of target photons in #/cm3/Hz
    double Phi;        // spectrum of electrons/positrons
    double Ep;         // proton energy Ep=x/Ee
    double fph_jet;    // number density of target photons in #/cm3/Hz

    Ep = E / pow(10., x);
    fp = colliding_protons(spline_Jp, acc_Jp, gp_min, gp_max, Ep);
    fph_jet = photons_jet(eta, Ep, spline_ng, acc_ng, nu_min, nu_max);
    fph = fph_jet;    // all the target photons are now included into the
                      // fph_jet function
    //	double Tstar = 2.7;
    //	fph = photon_field(eta,Ep,Tstar);
    Phi = PhiFunc(eta, eta_zero, pow(10., x), product);
    return fp * fph * Phi * log(10.) / Ep;
}

//************************************************************************************************************
double PhiFunc(double eta, double eta0, double x, std::string product) {
    // eqs 27,28,29 etc for gamma rays with interpolation in the tables given by
    // KA08 and eqs 31,32,33,34,35,36,37,38,39,40,41 and tables for leptons

    double Phi;                  // spectrum of products
    double xminus, xplus;        // eq. 19 from KA08 for min/max energy of pion
    double r = 0.146;            // mpion/mproton = 0.146 for above expressions
    double s, delta, Beta;       // the parameters for spectrum
    double xeta = eta / eta0;    // x_eta to distinguish from x = Ep/Ei -- it's
                                 // ρ sometimes in KA08

    tables_photomeson(s, delta, Beta, product, xeta);

    xplus =
        1. / (2. + 2. * eta) *
        (eta + r * r + sqrt((eta - r * r - 2. * r) * (eta - r * r + 2. * r)));
    xminus =
        1. / (2. + 2. * eta) *
        (eta + r * r - sqrt((eta - r * r - 2. * r) * (eta - r * r + 2. * r)));

    if (product.compare("gamma_rays") == 0) {
        double y;
        y = (x - xminus) / (xplus - xminus);
        if (x > xminus && x < xplus) {
            Phi = Beta * exp(-s * pow(log(x / xminus), delta)) *
                  pow(log(2. / (1. + y * y)), 2.5 + 0.4 * log(eta / eta0));
        } else if (x <= xminus) {
            Phi = Beta * pow(log(2.), 2.5 + 0.4 * log(eta / eta0));
        } else {
            Phi = 1.e-100;
        }
    } else if (product.compare("electrons") == 0 ||
               product.compare("antielectron") == 0) {
        double yprime, xplusprime, xminusprime, psi;
        xplus = 1. / (2. + 2. * eta) *
                (eta - 2. * r +
                 sqrt(eta * (eta - 4. * r * (1. + r))));    // eq40 max
        xminus = 1. / (2. + 2. * eta) *
                 (eta - 2. * r -
                  sqrt(eta * (eta - 4. * r * (1. + r))));    // eq40 min
        xplusprime = xplus;                                  // eq40 comments
        xminusprime = xminus / 2.;                           // eq40 comments
        if (xeta >= 4.0) {
            psi = 6. * (1. - exp(1.5 * (4. - xeta)));
        } else {
            psi = 0.;
        }

        if (x < xminusprime) {    // eq33
            Phi = Beta * pow(log(2.), psi);
        } else if (xminusprime <= x and x < xplusprime) {    // eq31
            yprime = (x - xminusprime) / (xplusprime - xminusprime);
            Phi = Beta * exp(-s * pow(log(x / xminusprime), delta)) *
                  pow(log(2. / (1. + yprime * yprime)), psi);
        } else {
            Phi = 1.e-100;
        }
    } else if (product.compare("positrons") == 0 ||
               product.compare("antimuon") == 0 ||
               product.compare("electron") == 0) {
        double yprime, xplusprime, xminusprime, psi;
        xplusprime = xplus;
        xminusprime = xminus / 4.;

        psi = 2.5 + 1.4 * log(xeta);

        if (x < xminusprime) {    // eq33
            Phi = Beta * pow(log(2.), psi);
        } else if (xminusprime <= x and x < xplusprime) {    // eq31
            yprime = (x - xminusprime) / (xplusprime - xminusprime);
            Phi = Beta * exp(-s * pow(log(x / xminusprime), delta)) *
                  pow(log(2. / (1. + yprime * yprime)), psi);
        } else {
            Phi = 1.e-200;
        }
    } else if (product.compare("muon") == 0) {
        double yprime, xplusprime, xminusprime, psi;
        if (xeta < 2.14) {
            xplusprime = 0.427 * xplus;
        } else if (xeta >= 2.14 and xeta < 10.) {
            xplusprime = (0.427 + 0.0729 * (xeta - 2.14)) * xplus;
        } else {
            xplusprime = xplus;
        }
        xminusprime = 0.427 * xminus;

        psi = 2.5 + 1.4 * log(xeta);

        if (x < xminusprime) {    // eq33
            Phi = Beta * pow(log(2.), psi);
        } else if (xminusprime <= x and x < xplusprime) {    // eq31
            yprime = (x - xminusprime) / (xplusprime - xminusprime);
            Phi = Beta * exp(-s * pow(log(x / xminusprime), delta)) *
                  pow(log(2. / (1. + yprime * yprime)), psi);
        } else {
            Phi = 1.e-200;
        }
    } else {
        Phi = 1.e-100;
    }

    return Phi;
}
// The tables from KA08 for photomeson that give s,δ and B
void tables_photomeson(double &s, double &delta, double &Beta,
                       std::string product, double xeta) {
    // Gamma rays from neutral pion decay:
    int sizeTable;
    if (product.compare("electrons") == 0 ||
        product.compare("antielectron") == 0) {
        sizeTable = 10;
    } else if (product.compare("positrons") == 0 ||
               product.compare("muon") == 0 ||
               product.compare("antimuon") == 0 ||
               product.compare("electron") == 0) {
        sizeTable = 20;
    } else {
        std::cout << "wrong population" << std::endl;
        sizeTable = 2;
    }
    gsl_interp_accel *acc_sigma = gsl_interp_accel_alloc();
    gsl_spline *spline_sigma = gsl_spline_alloc(gsl_interp_cspline, sizeTable);
    gsl_interp_accel *acc_delta = gsl_interp_accel_alloc();
    gsl_spline *spline_delta = gsl_spline_alloc(gsl_interp_cspline, sizeTable);
    gsl_interp_accel *acc_Beta = gsl_interp_accel_alloc();
    gsl_spline *spline_Beta = gsl_spline_alloc(gsl_interp_cspline, sizeTable);

    if (product.compare("electrons") == 0) {
        gsl_spline_init(spline_sigma, etaeTable, seTable, sizeTable);
        gsl_spline_init(spline_delta, etaeTable, deltaeTable, sizeTable);
        gsl_spline_init(spline_Beta, etaeTable, BetaeTable, sizeTable);
    }    // positrons from charged pion decay:
    else if (product.compare("positrons") == 0) {
        gsl_spline_init(spline_sigma, etaposTable, sposTable, sizeTable);
        gsl_spline_init(spline_delta, etaposTable, deltaposTable, sizeTable);
        gsl_spline_init(spline_Beta, etaposTable, BetaposTable, sizeTable);
    }    // muon neutrinos
    else if (product.compare("muon") == 0) {
        gsl_spline_init(spline_sigma, eta_muonTable, s_muonTable, sizeTable);
        gsl_spline_init(spline_delta, eta_muonTable, delta_muonTable,
                        sizeTable);
        gsl_spline_init(spline_Beta, eta_muonTable, Beta_muonTable, sizeTable);
    }    // muon antineutrinos
    else if (product.compare("antimuon") == 0) {
        gsl_spline_init(spline_sigma, eta_antimuonTable, s_antimuonTable,
                        sizeTable);
        gsl_spline_init(spline_delta, eta_antimuonTable, delta_antimuonTable,
                        sizeTable);
        gsl_spline_init(spline_Beta, eta_antimuonTable, Beta_atnimuonTable,
                        sizeTable);
    }    // electron neutrinos
    else if (product.compare("electron") == 0) {
        gsl_spline_init(spline_sigma, eta_electronTable, s_electronTable,
                        sizeTable);
        gsl_spline_init(spline_delta, eta_electronTable, delta_electronTable,
                        sizeTable);
        gsl_spline_init(spline_Beta, eta_electronTable, Beta_electronTable,
                        sizeTable);
    }    // electron antineutrinos
    else if (product.compare("antielectron") == 0) {
        gsl_spline_init(spline_sigma, eta_antielectronTable,
                        s_antielectronTable, sizeTable);
        gsl_spline_init(spline_delta, eta_antielectronTable,
                        delta_antielectronTable, sizeTable);
        gsl_spline_init(spline_Beta, eta_antielectronTable,
                        Beta_antielectronTable, sizeTable);
    }

    s = gsl_spline_eval(spline_sigma, xeta, acc_sigma);
    delta = gsl_spline_eval(spline_delta, xeta, acc_delta);
    Beta = gsl_spline_eval(spline_Beta, xeta, acc_Beta);

    gsl_spline_free(spline_sigma), gsl_interp_accel_free(acc_sigma);
    gsl_spline_free(spline_delta), gsl_interp_accel_free(acc_delta);
    gsl_spline_free(spline_Beta), gsl_interp_accel_free(acc_Beta);
}

}    // namespace kariba
