import sys
sys.path.insert(0,"./src/wrapper/")
import py_Delta_Sigma_miscentering as pyDSm
import matplotlib.pyplot as plt
import numpy as np
plt.rc("font", size=12)


R, sigma = np.genfromtxt("test_data/test.txt", unpack=True)
cosmo = {"h":0.7,"om":0.3,"ok":0.0}
cosmo["ode"]=1.0-cosmo["om"]

params = {"Mass": 3*10**14,"NR":300,"Rmin":0.01,
                "Rmax":200.0,"Nbins":15,"R_bin_min":0.01,"R_bin_max":200.0,
                "delta":200,"Rmis":0.25,"fmis":0.25,
                "averaging":1,"single":1}
params["concentration"] = 5.0 #Arbitrary

result = pyDSm.calc_Delta_Sigma_miscentering(R, sigma, cosmo, params)

import pickle
R_old, sigma_old, sigma_miscentered_old, sigma_m_old,\
    miscentered_delta_sigma_old, delta_sigma_single_old,\
    ave_miscentered_delta_sigma_old,\
    ave_delta_sigma_single_old = pickle.load(open("test_data/output.p", "r"))

R = result["R"]
sigma = result['sigma']
sigma_miscentered = result['miscentered_sigma']
sigma_m = result['sigma_single']
miscentered_delta_sigma = result['miscentered_delta_sigma']
delta_sigma_single = result['delta_sigma_single']
ave_miscentered_delta_sigma = result['ave_miscentered_delta_sigma']
ave_delta_sigma_single = result['ave_delta_sigma_single']
np.testing.assert_array_equal(R, R_old)
np.testing.assert_array_equal(sigma, sigma_old)
np.testing.assert_array_equal(sigma_miscentered, sigma_miscentered_old)
np.testing.assert_array_equal(sigma_m, sigma_m_old)
np.testing.assert_array_equal(miscentered_delta_sigma, miscentered_delta_sigma_old)
np.testing.assert_array_equal(delta_sigma_single, delta_sigma_single_old)
np.testing.assert_array_equal(ave_miscentered_delta_sigma, ave_miscentered_delta_sigma_old)
np.testing.assert_array_equal(ave_delta_sigma_single, ave_delta_sigma_single_old)

plt.loglog(R,sigma,label=r"$\Sigma$")
plt.loglog(R,sigma_miscentered,label=r"$\Sigma_{\rm mis}$",ls='--')
plt.loglog(R,sigma_m,label=r"$\Sigma(R|R_{\rm mis})$")
plt.legend()
plt.xlabel(r"$R\ [{\rm Mpc}/h]$",fontsize=24)
plt.ylabel(r"$\Sigma\ [{\rm M_\odot}\ h/{\rm pc^2}]$",fontsize=24)
plt.subplots_adjust(bottom=0.15, left=0.15)
plt.show()
plt.clf()
