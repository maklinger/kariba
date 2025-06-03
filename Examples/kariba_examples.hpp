#include <stdarg.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <fenv.h>
#include <omp.h>

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>


using namespace std;

void plot_write(int size,double *en,double *lum,char path[],double dist,double redshift);					
void plot_write(int size,const double *en,const double *lum,char path[],double dist,double redshift);	
void plot_write(int size,const double *p,const double *g,const double *pdens,const double *gdens,
				char path[]);

void sum_zones(int size_in,int size_out,double* input_en,double* input_lum,double* en,double* lum);
void sum_ext(int size_in,int size_out,const double* input_en,const double* input_lum,double* en,double* lum);
double integrate_lum(int size,double numin,double numax,const double* input_en,const double* input_lum);
double photon_index(int size,double numin,double numax,const double* input_en,const double* input_lum);

void clean_file(char path[],bool check);
void read_params(string file,double *pars);





