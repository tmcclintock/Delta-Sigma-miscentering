/*
  This is an interface between the DeltaSigma code and 
  a piece of code that can pass in a cosmology and 
  a power spectrum.
*/
#include "../sigma_mis/sigma_mis.h"
#include "../delta_sigma_mis/delta_sigma_mis.h"
#include "../miscentered_sigma/miscentered_sigma.h"
#include "../miscentered_delta_sigma/miscentered_delta_sigma.h"
#include "../ave_miscentered_delta_sigma/ave_miscentered_delta_sigma.h"
#include "../constants/constants.h"
#include "../cosmology/cosmology.h"

#ifndef INTERFACE
#define INTERFACE
typedef struct interface_parameters{
  double Mass;
  double concentration;
  double Rmis;
  int delta;
  int miscentering;
  int averaging;
  int single;
  int Nbins;
  double R_bin_min;
  double R_bin_max;
}interface_parameters;
#endif

#ifndef WRAPPER_OUTPUT
#define WRAPPER_OUTPUT
typedef struct wrapper_output{
  double*R;
  double*sigma;
  double*Rbins;
  double*sigma_mis;
  double*delta_sigma_mis;
  double*miscentered_sigma;
  double*miscentered_delta_sigma;
  double*ave_miscentered_delta_sigma;
  double*ave_delta_sigma_mis;
}wrapper_output;
#endif

int interface(int NR, cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs);

int python_interface(int NR, double h, double om,
		     double Mass, double concentration,
		     double Rmis, int delta,
		     int averaging, int single,
		     int Nbins,
		     double R_bin_min, double R_bin_max,
		     double*R,
		     double*sigma,
		     double*Rbins,
		     double*sigma_mis,
		     double*delta_sigma_mis,
		     double*miscentered_sigma,
		     double*miscentered_delta_sigma,
		     double*ave_miscentered_delta_sigma,
		     double*ave_delta_sigma_mis);

