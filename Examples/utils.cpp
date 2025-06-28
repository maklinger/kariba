#include <kariba/Radiation.hpp>
#include <kariba/constants.hpp>

#include "kariba_examples.hpp"

namespace karcst = kariba::constants;

// Read input parameters from file, store them in an array
void read_params(std::string file, double *pars) {
    std::ifstream inFile;
    inFile.open(file.c_str());
    std::string line;
    int line_nb = 0;
    if (!inFile) {
        std::cerr << "Can't open input file\n";
        exit(1);
    }
    while (getline(inFile, line)) {
        // Remove whitespace from the beginning of the line
        line.erase(line.begin(),
                   std::find_if(line.begin(), line.end(), [](unsigned char c) {
                       return !std::isspace(c);
                   }));
        if (line[0] == '#') {
            continue;
        } else {
            pars[line_nb] = atof(line.c_str());
            line_nb++;
        }
    }
    inFile.close();
}

// Plots a given array in units of ergs (x axis) and erg/s/Hz (y axis) to the
// file on the provided path. Note: the factor 1+z in
// the specific luminosity calculation is to ensure that the output spectrum
// only moves to lower frequency, not up/down.
void plot_write(int size, const double *en, const double *lum,
                const std::string &path, double redsh) {
    std::ofstream file;
    file.open(path, std::ios::app);

    for (int k = 0; k < size; k++) {
        file << en[k] / (karcst::herg * (1. + redsh)) << " " << lum[k] << "\n";
    }

    file.close();
}

void plot_write(int size, const double *p, const double *g, const double *pdens,
                const double *gdens, const std::string &path) {

    std::ofstream file;
    file.open(path, std::ios::app);
    for (int k = 0; k < size; k++) {
        file << p[k] << " " << g[k] << " " << pdens[k] << " " << gdens[k]
             << "\n";
    }

    file.close();
}

// Used for summing individual zone contributions for a generic spectral
// component from code: pre/post particle acceleration synchrotron, pre/post
// particle acceleration Comptonization The second function does the same, but
// sums the disk/corona/bb to the total jet spectrum. The reason for the const
// arryas in input is that the input arrays are directly accessed from the
// ShSDisk class, which are const
void sum_zones(int size_in, int size_out, double *input_en, double *input_lum,
               double *en, double *lum) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_akima, size_in);
    gsl_spline_init(input_spline, input_en, input_lum, size_in);

    for (int i = 0; i < size_out; i++) {
        if (en[i] > input_en[0] && en[i] < input_en[size_in - 1]) {
            lum[i] = lum[i] + gsl_spline_eval(input_spline, en[i], acc);
        }
    }
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

void sum_ext(int size_in, int size_out, const double *input_en,
             const double *input_lum, double *en, double *lum) {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *input_spline = gsl_spline_alloc(gsl_interp_akima, size_in);
    gsl_spline_init(input_spline, input_en, input_lum, size_in);

    for (int i = 0; i < size_out; i++) {
        if (en[i] > input_en[0] && en[i] < input_en[size_in - 1]) {
            lum[i] = lum[i] + gsl_spline_eval(input_spline, en[i], acc);
        }
    }
    gsl_spline_free(input_spline), gsl_interp_accel_free(acc);
}

// Simple numerical integral to calculate the luminosity between numin and numax
// of a given array; input units are erg for the frequency/energy array,
// erg/s/Hz for the luminosity array to be integrated, Hz for the integration
// bounds; size is the dimension of the input arrays Note: this uses a VERY
// rough method and wide bins, so thread carefully
double integrate_lum(int size, double numin, double numax,
                     const double *input_en, const double *input_lum) {
    double temp = 0.;
    for (int i = 0; i < size - 1; i++) {
        if (input_en[i] / karcst::herg > numin &&
            input_en[i + 1] / karcst::herg < numax) {
            temp = temp + (1. / 2.) *
                              (input_en[i + 1] / karcst::herg -
                               input_en[i] / karcst::herg) *
                              (input_lum[i + 1] + input_lum[i]);
        }
    }
    return temp;
}

// Overly simplified estimate of the photon index between numin and numax of a
// given array; input is the same as integrate_lum. Note that this assumes
// input_lum is a power-law in shape
double photon_index(int size, double numin, double numax,
                    const double *input_en, const double *input_lum) {
    int counter_1 = 0, counter_2 = 0;
    double delta_y, delta_x, gamma;
    for (int i = 0; i < size; i++) {
        if (input_en[i] / karcst::herg < numin) {
            counter_1 = i;
        }
        if (input_en[i] / karcst::herg < numax) {
            counter_2 = i;
        }
    }
    delta_y = log10(input_lum[counter_2]) - log10(input_lum[counter_1]);
    delta_x = log10(input_en[counter_2] / karcst::herg) -
              log10(input_en[counter_1] / karcst::herg);
    gamma = delta_y / delta_x - 1.;
    return gamma;
}

// Prepares files for above printing functions at the start of the run. There
// are two reasons this exists: 1) specifying the units of the output obviously
// makes things easier to read and 2) S-lang does not allow to pass ofstream
// objects to functions, so we need to pass a path and open the file from the
// path inside the write function. As a result, it's impossible to just truncate
// and clean the files at the start of each iteration. This is only relevant for
// the cyclosyn_zones, compton_zones, and numdens files
void clean_file(const std::string &path, bool check) {
    std::ofstream file;
    file.open(path, std::ios::trunc);

    if (check == true) {
        file << std::left << std::setw(20) << "#nu [Hz] " << "Flux [mJy]\n";
    } else {
        file << std::left << std::setw(20) << "#p [g cm s-1] " << std::setw(20)
             << "g [] " << std::setw(20) << " n(p) [# cm^-3 p^-1] "
             << " n(g) [# cm^-3 g^-1]\n";
    }

    file.close();
}
