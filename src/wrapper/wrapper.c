#include "wrapper.h"

int interface(int NR, cosmology cosmo, interface_parameters*params,
	      wrapper_output*outputs){

  int i;

  double*R=outputs->R;
  double*sigma=outputs->sigma;
  double*Rbins=outputs->Rbins;
  double*sigma_mis=outputs->sigma_mis;
  double*delta_sigma_mis=outputs->delta_sigma_mis;
  double*miscentered_sigma=outputs->miscentered_sigma;
  double*miscentered_delta_sigma=outputs->miscentered_delta_sigma;
  double*ave_miscentered_delta_sigma=outputs->ave_miscentered_delta_sigma;
  double*ave_delta_sigma_mis=outputs->ave_delta_sigma_mis;

  //Used to hold integration errors for some routines
  double*err=(double*)malloc(NR*sizeof(double));

  double Mass=params->Mass;
  double concentration=params->concentration;
  double Rmis=params->Rmis;
  int delta=params->delta;

  int Nbins=params->Nbins;
  double R_bin_min=params->R_bin_min;
  double R_bin_max=params->R_bin_max;

  int averaging=params->averaging;
  int single=params->single;

  calc_miscentered_sigma(R,Mass,concentration,delta,Rmis,R,sigma,NR,miscentered_sigma,err,cosmo);
  calc_miscentered_delta_sigma(R,Mass,concentration,delta,Rmis,R,sigma,miscentered_sigma,NR,miscentered_delta_sigma,err,cosmo);

  if(single){
    calc_sigma_mis(R,Mass,concentration,delta,Rmis,R,sigma,NR,sigma_mis,err,cosmo);
    calc_delta_sigma_mis(R,Mass,concentration,delta,Rmis,R,sigma,sigma_mis,NR,delta_sigma_mis,err,cosmo);
    if(averaging){
      calc_ave_miscentered_delta_sigma(R,NR,delta_sigma_mis,Nbins,R_bin_min,R_bin_max,Rbins,ave_delta_sigma_mis);
    }
  }

  if(averaging){
    calc_ave_miscentered_delta_sigma(R,NR,miscentered_delta_sigma,Nbins,R_bin_min,R_bin_max,Rbins,ave_miscentered_delta_sigma);
  }

  free(err);
  return 0;
}

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
		     double*ave_delta_sigma_mis){
  cosmology*cosmo = (cosmology*)malloc(sizeof(cosmology));
  cosmo->H0 = h*100.;
  cosmo->h = h;
  cosmo->om = om;

  interface_parameters*params=
    (interface_parameters*)malloc(sizeof(interface_parameters));
  params->Mass=Mass;
  params->concentration=concentration;
  params->delta=delta;
  params->Rmis=Rmis;
  params->averaging=averaging; //1 is true
  params->single=single; //1 is true
  params->Nbins=Nbins;
  params->R_bin_min=R_bin_min;
  params->R_bin_max=R_bin_max;

  wrapper_output*outputs=(wrapper_output*)malloc(sizeof(wrapper_output));
  outputs->R=R;
  outputs->sigma=sigma;
  outputs->sigma_mis=sigma_mis;
  outputs->delta_sigma_mis=delta_sigma_mis;
  outputs->miscentered_sigma=miscentered_sigma;
  outputs->miscentered_delta_sigma=miscentered_delta_sigma;
  outputs->Rbins=Rbins;
  outputs->ave_miscentered_delta_sigma=ave_miscentered_delta_sigma;
  outputs->ave_delta_sigma_mis=ave_delta_sigma_mis;

  interface(NR,*cosmo,params,outputs);

  free(cosmo);
  free(params);
  free(outputs);
  return 0;
}
		     
