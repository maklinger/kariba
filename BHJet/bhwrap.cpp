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

void read_params(std::string file, double *pars);

extern void jetinterp(double *ear, double *energ, double *phot, double *photar,
                      int ne, int newne);
extern void jetmain(double *ear, int ne, double *param, double *photeng,
                    double *photspec);

int main() {

    double start = omp_get_wtime();

    int npar = 28;
    int ne = 201;
    double emin = -10;
    double emax = 10;
    double einc = (emax - emin) / ne;

    double *ebins = new double[ne]();
    double *param = new double[npar]();
    double *spec = new double[ne - 1]();
    double *dumarr = new double[ne - 1]();

    for (int i = 0; i < ne; i++) {
        ebins[i] = pow(10, (emin + i * einc));
    }

    read_params("Input/ip.dat", param);

    jetmain(ebins, ne - 1, param, spec, dumarr);

    delete[] ebins, delete[] param, delete[] spec, delete[] dumarr;

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
