"""
Created on Tue Apr 20th, 2021

Auswertung des N-body Problems
@author: Martin, Florian, Sebastian
"""

import numpy as np
import matplotlib.pyplot as plt

Daten_rk4 = np.loadtxt("rk4-solution.csv",delimiter=';')
# Daten_fwd = np.loadtxt("fwd-solution.csv")
# Daten_lf = np.loadtxt("lf-solution.csv")

plt.figure(dpi=300)
plt.plot(Daten_rk4[:,0], Daten_rk4[:,1], label='Spektrum')
plt.legend()
# plt.savefig("0421_Integration_rk4.pdf")