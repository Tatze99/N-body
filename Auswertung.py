"""
Created on Tue Apr 20th, 2021

Auswertung des N-body Problems
@author: Martin, Florian, Sebastian
"""

import numpy as np
import matplotlib.pyplot as plt

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

command = "rk4"
Daten = np.loadtxt(command+"-solution.csv",delimiter=';')
mass = np.loadtxt("Input.csv",delimiter=';',usecols=[6])
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto']

steps = len(Daten[:,0])
number = 10
n = 10
m = mass[0:n]
names = Namen[0:n]
t  = Daten[:,0]
x  = Daten[:,0*n+1:1*n+1]
y  = Daten[:,1*n+1:2*n+1]
z  = Daten[:,2*n+1:3*n+1]
vx = Daten[:,3*n+1:4*n+1]
vy = Daten[:,4*n+1:5*n+1]
vz = Daten[:,5*n+1:6*n+1]

#%%
# Plot the trajectories
plt.figure(dpi=300)
plt.plot(x[:,0:number], y[:,0:number], label=names[0:number])
plt.legend()
plt.savefig("0421_Integration_"+command+".pdf")

#%%
# Plot the energy
kin = kinetic_energy(vx, vy, vz, m)
pot = potential_energy(x, y, z, m)

plt.figure(dpi=300)
plt.plot(t, kin, label='kinetic energy')
plt.plot(t, pot, label='potential energy')
plt.plot(t, kin + pot, label='total energy')
plt.legend()
