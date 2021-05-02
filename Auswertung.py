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
time  = Daten[:,0]
number = 5          # number of planets to display
n = len(mass)-1     # total number of planets


# create the variables and assign them their values via a loop
var_names = ["x", "y", "z","vx", "vy", "vz"]
for i,name in enumerate(var_names):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]

#%%
# Plot the trajectories
plt.figure(dpi=300)
plt.plot(x[:,0:number], y[:,0:number],'-',linewidth=1, label=Namen[0:number])
plt.legend()
plt.savefig("0421_Integration_"+command+".pdf")

#%%
# Plot the energy
kin = kinetic_energy(vx, vy, vz, mass)
pot = potential_energy(x, y, z, mass)
# kin = Daten[:,6*n+1]
# pot = Daten[:,6*n+2]

plt.figure(dpi=300)
legend = ['kinetic energy', 'potential energy', 'total energy']
for i,y in enumerate([kin, pot, kin+pot]):
    plt.plot(time,y,label=legend[i])
plt.legend()