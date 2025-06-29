#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fenv.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

void read_params(std::string file, std::vector<double> &pars);

extern void jetinterp(std::vector<double> &ear, std::vector<double> &energ,
                      std::vector<double> &phot, std::vector<double> &photar,
                      size_t ne, size_t newne);
extern void jetmain(std::vector<double> &ear, size_t ne,
                    std::vector<double> &param, std::vector<double> &photeng,
                    std::vector<double> &photspec);

int main() {

    double start = omp_get_wtime();

    size_t npar = 28;
    size_t ne = 201;
    double emin = -10;
    double emax = 10;
    double einc = (emax - emin) / ne;

    std::vector<double> ebins(ne, 0.0);
    std::vector<double> param(npar, 0.0);
    std::vector<double> spec(ne - 1, 0.0);
    std::vector<double> dumarr(ne - 1, 0.0);

    for (size_t i = 0; i < ne; i++) {
        ebins[i] = pow(10, (emin + i * einc));
    }

    read_params("Input/ip.dat", param);

    jetmain(ebins, ne - 1, param, spec, dumarr);

    double end = omp_get_wtime();
    std::cout << "Total running time: " << end - start << " seconds\n";

    // system("python3 Plot_separate.py");

    return EXIT_SUCCESS;
}    // ----------  end of function main  ----------

// Read input file
//
//  @param file		Parameters file
//
//  @return pars         Parameters
//
void read_params(std::string file, std::vector<double> &pars) {
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
