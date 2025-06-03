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

using namespace std;

extern void jetmain(double *ear, int ne, double *param, double *photeng,
                    double *photspec);

extern "C" void pyjetmain(double *ear, int ne, double *param, double *photeng,
                          double *photspec) {
    jetmain(ear, ne, param, photeng, photspec);
} // ----------  end of function main  ----------
