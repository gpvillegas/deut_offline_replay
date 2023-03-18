import numpy as np
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
import matplotlib.pyplot as plt

# runs: 21049 - 21058
I_avg_LD2     = np.array([9.965, 14.828, 24.823, 34.051, 38.137, 41.540, 42.459, 51.005, 64.827, 67.364 ]) # this is really yield, I just could not use the same name twice
LD2_counts = np.array([276149.0, 276019.0, 565177.0, 359402.0, 394579.0, 434953.0, 742710.0, 411562.0, 467883.0, 559969.])
LD2_Q = np.array([5.561, 5.577, 12.100, 14.693, 16.217, 17.940, 30.685, 25.560, 29.575, 36.718])
LD2_etrk_eff = np.array([0.987, 0.986, 0.985, 0.984, 0.983, 0.983, 0.983, 0.981, 0.981, 0.979])
LD2_tLT = np.array([1.000, 0.998, 0.953, 1.006, 1.006, 1.006, 1.005, 1.014, 0.986, 0.980])
PS_factor = np.array([1., 1., 1., 2., 2., 2., 2., 3., 3., 3.])

LD2_Y = (LD2_counts*PS_factor) / (LD2_Q*LD2_etrk_eff*LD2_tLT)
LD2_Y_err = np.sqrt(LD2_counts) *PS_factor / (LD2_Q*LD2_etrk_eff*LD2_tLT)


LD2_rel_yield = LD2_Y/LD2_Y[0]
LD2_rel_yield_err = LD2_Y_err/LD2_Y[0]


#LD2_Y = np.array([50317.504, 50294.898, 49741.668, 49431.715, 49167.832, 49059.602, 49012.387, 48545.090, 49105.715, 47643.395 ])

#LD2_Y_err = np.sqrt(LD2_Y)
#LD2_yield = unumpy.uarray(LD2_counts, LD2_counts_err)

# runs 21059, 21060, 21061
I_avg_C12 = np.array([14.858, 29.442, 39.796, 55.328, 64.626])
C12_counts = np.array([162756.0, 266221., 286626., 613042.0, 600152.0])
C12_Q = np.array([9.835, 16.009, 17.258, 37.007, 36.588   ])
C12_etrk_eff = np.array([0.987, 0.988, 0.987, 0.986, 0.986])
C12_tLT = np.array([1.000, 1.001, 1.000, 1.0, 0.996])
PS_factor_c12 = np.array([1,1,1,1,1])

C12_Y = (C12_counts*PS_factor_c12) / (C12_Q*C12_etrk_eff*C12_tLT)
C12_Y_err = np.sqrt(C12_counts) *PS_factor_c12 / (C12_Q*C12_etrk_eff*C12_tLT)

C12_rel_yield = C12_Y/C12_Y[0]
C12_rel_yield_err = C12_Y_err/C12_Y[0]

#------------




plt.errorbar(I_avg_C12, C12_rel_yield, C12_rel_yield_err, marker='o', markersize=10, linestyle='None',  color='k', label='C12')
plt.errorbar(I_avg_LD2, LD2_rel_yield, LD2_rel_yield_err, marker='^', markersize=10, linestyle='None', color='r', label='LD2' )

plt.title('Preliminary d(e,e\'p) Luminosity Scans (March 18, 2023)',  fontsize=18)
plt.xlabel('Average BCM4A Current [uA]', fontsize=18)
plt.ylabel('Relative Charge-Normalized Yield', fontsize=18)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)

plt.grid()

plt.legend(fontsize=20)

#plt.errorbar(I_avg_LD2, LD2_counts, LD2_counts_err, marker='o', color='k' )

#plt.plot(I_avg_LD2, LD2_T2_scl_counts_norm, marker='o', color='k' )




plt.show()
