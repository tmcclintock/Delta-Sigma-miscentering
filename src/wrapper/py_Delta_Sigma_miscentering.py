"""
This is the python wrapper that calls the python_wrapper()
c function. This interfaces through c_types so that the user
doesn't have to.
"""
import numpy as np
from ctypes import c_double,c_int,POINTER,cdll
import inspect
import os
library_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))+"/Delta_Sigma_miscentering.so"
dslib = cdll.LoadLibrary(library_path)
interface = dslib.python_interface
interface.restype = c_int
"""
Arguments to the interface are: 
NR,h,om,
Mass,concentration,
Rmis,delta,
single,
averaging, Nbins,
R_bin_min,R_bin_max,
R,
sigma,
Rbins,
sigma_single,
delta_sigma_single
miscentered_sigma,
miscentered_delta_sigma,
ave_miscentered_delta_sigma
ave_delta_sigma_single
"""
interface.argtypes=[c_int, c_double, c_double,
                    c_double, c_double,
                    c_double, c_int,
                    c_int,
                    c_int,c_int,
                    c_double,c_double,
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double),
                    POINTER(c_double)]

def calc_Delta_Sigma_miscentering(R, sigma, cosmo_dict, params):
    """Calculates the DeltaSigma profile given some cosmology, matter power spectra, and input parameters (e.g. mass, concentraton, etc.)

    Note: Mass units are Msun/h. Distances are Mpc/h comoving.

    R (array_like): Radii of input surface mass density; Mpc/h
    sigma (array_like): Surface mass density; Msun h/pc^2
    cosmo_dict (dictionary): Contains key-value pairs of cosmological parameters. Required parameters: h, om, and ode.
    params (dictionary): Contains key-value pairs of halo parameters, including: Mass, delta, Rmis, fmis, concentration, NR, Nbins, R_bin_min, R_bin_max, averaging, single (for a single miscentered halo).

    Returns:
        output (dictionary): Contains key-value pairs of all possible quantities assosciated with the halo.
    """
    Mass,concentration,delta = params["Mass"],params["concentration"],params['delta']
    NR = params["NR"]
    Rmis,fmis = params["Rmis"], params["fmis"]

    if "single" in params: single = params['single']
    else: single = 0
    if "averaging" in params: 
        averaging = params['averaging']
        Nbins = params['Nbins']
        R_bin_min = params['R_bin_min']
        R_bin_max = params['R_bin_max']
    else: #Default values to pass to C
        averaging = 0
        Nbins = 2
        R_bin_min = min(R)
        R_bin_max = max(R)

    h,om = cosmo_dict['h'],cosmo_dict['om']

    R = R.astype("float64")
    sigma = sigma.astype("float64")
    R_in = R.ctypes.data_as(POINTER(c_double))
    sigma_in = sigma.ctypes.data_as(POINTER(c_double))
    sigma_single = np.zeros(NR)
    sigma_single_in = sigma_single.ctypes.data_as(POINTER(c_double))
    delta_sigma_single = np.zeros(NR)
    delta_sigma_single_in = delta_sigma_single.ctypes.data_as(POINTER(c_double))
    miscentered_sigma = np.zeros(NR)
    miscentered_sigma_in = miscentered_sigma.ctypes.data_as(POINTER(c_double))
    miscentered_delta_sigma = np.zeros(NR)
    miscentered_delta_sigma_in = miscentered_delta_sigma.ctypes.data_as(POINTER(c_double))

    Rbins = np.zeros(Nbins)
    Rbins_in = Rbins.ctypes.data_as(POINTER(c_double))
    ave_miscentered_delta_sigma = np.zeros(Nbins)
    ave_miscentered_delta_sigma_in = ave_miscentered_delta_sigma.ctypes.data_as(POINTER(c_double))
    ave_delta_sigma_single = np.zeros(Nbins)
    ave_delta_sigma_single_in = ave_delta_sigma_single.ctypes.data_as(POINTER(c_double))

    result = interface(NR,h,om,
                       Mass,concentration,
                       Rmis,delta,
                       averaging,single,
                       Nbins,R_bin_min,R_bin_max,
                       R_in,
                       sigma_in,
                       Rbins_in,
                       sigma_single_in,
                       delta_sigma_single_in,
                       miscentered_sigma_in,
                       miscentered_delta_sigma_in,
                       ave_miscentered_delta_sigma_in,
                       ave_delta_sigma_single_in)

    #Now build a dictionary and return it
    output = {"R":R,"sigma":sigma}
    output["miscentered_sigma"] = miscentered_sigma
    output["miscentered_delta_sigma"] = miscentered_delta_sigma
    
    if single:
        output["sigma_single"] = sigma_single
        output["delta_sigma_single"] = delta_sigma_single
    
    if averaging:
        output["Rbins"] = Rbins
        output["ave_miscentered_delta_sigma"] = ave_miscentered_delta_sigma
        if single:
            output["ave_delta_sigma_single"] = ave_delta_sigma_single

    return output
