"""
Created on Tue Apr 20th, 2021

Auswertung des N-body Problems
@author: Martin, Florian, Sebastian
"""

import numpy as np
import matplotlib.pyplot as plt

# Initial = np.loadtxt("Initial.csv", delimiter=';')
Daten_rk4 = np.loadtxt("rk4-solution.csv",delimiter=';')
# Daten_fwd = np.loadtxt("fwd-solution.csv")
# Daten_lf = np.loadtxt("lf-solution.csv")

plt.figure(dpi=300)
plt.plot(Daten_rk4[:,1], Daten_rk4[:,3], label='Sonne')
plt.plot(Daten_rk4[:,2], Daten_rk4[:,4], label='Erde')
# plt.plot(Daten_rk4[:,0], Daten_rk4[:,1], label='Spektrum')
plt.legend()
# plt.savefig("0421_Integration_rk4.pdf")

steps = len(Daten_rk4[:,0])
m = [1,0.000003003]
n = 2

t = Daten_rk4[:,0]
x = Daten_rk4[:,0*n+1:1*n+1]
y = Daten_rk4[:,1*n+1:2*n+1]
z = Daten_rk4[:,2*n+1:3*n+1]
vx = Daten_rk4[:,3*n+1:4*n+1]
vy = Daten_rk4[:,4*n+1:5*n+1]
vz = Daten_rk4[:,5*n+1:6*n+1]


def kinetic_energy(vx, vy, vz, m):
    energy = np.zeros(steps)
    for i in range(n):
        energy[:] += 0.5*m[i]*(vx[:,i]**2+vy[:,i]**2+vz[:,i]**2)
    return energy


def potential_energy(x, y, z, m):
    energy = np.zeros(steps)
    Matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            energy[:] += m[j]*m[i]/np.sqrt((x[:,i]-x[:,j])**2+(y[:,i]-y[:,j])**2+(z[:,i]-z[:,j])**2)
    return energy


kin = kinetic_energy(vx, vy, vz, m)
pot = potential_energy(x, y, z, m)

Energie = kin + pot

plt.figure(dpi=300)
plt.plot(t, kin, label='kinetic energy')
plt.plot(t, pot, label='potential energy')
plt.plot(t, Energie, label='total energy')
plt.legend()