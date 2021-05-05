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
number = 10          # number of planets to display
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
# plt.savefig("0421_Integration_"+command+".pdf")

#%%
# Plot the energy
kin = kinetic_energy(vx, vy, vz, mass)
pot = potential_energy(x, y, z, mass)
# kin = Daten[:,6*n+1]
# pot = Daten[:,6*n+2]

plt.figure(dpi=300)
legend = ['kinetic energy', 'potential energy', 'total energy']
for i,energy in enumerate([kin, pot, kin+pot]):
    plt.plot(time,energy,label=legend[i])
plt.legend()

#%%
# Plot in 3D
fig = plt.figure(figsize=(15, 6))
ax = plt.axes(projection='3d')

for i,name in enumerate(Namen):
    ax.plot(x[:,i], y[:,i], z[:,i], label=Namen[i])

ax.set_title('Trajectories of all planets')
ax.set_xlabel('$x$ in a.u.')
ax.set_ylabel('$y$ in a.u.')
ax.set_zlabel('$z$ in a.u.')
ax.legend()


#%%
from matplotlib import animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=120, metadata=dict(artist='Martin Beyer'), bitrate=-1)

fig = plt.figure(figsize=(6,3))
ax1 = plt.axes()
line, = ax1.plot([], [], lw=2)
plt.xlabel('$x$ in a.u.')
plt.ylabel('$y$ in a.u.')

ax1.plot(x[:,9],y[:,9],c='white')  # Cheap trick to get the plot dimensions right

lines = []
x2 = x[:,:number]
y2 = y[:,:number]
for index in range(number):
    lobj = ax1.plot(x[:1,index],y[:1,index],lw=2, label=Namen[index])[0]
    lines.append(lobj)
plt.legend(loc='right')
plt.xlim(-35,70)

def init():
    for line in lines:
        line.set_data([],[])
    return lines

def animate(i):
    for lnum,line in enumerate(lines):
        index = round(-0.012*lnum**6)
        if(i+index > 0):
            line.set_data(x[i+index:i,lnum], y[i+index:i,lnum]) # set data for each line separately.
        else:
            line.set_data(x[:i,lnum], y[:i,lnum])

    return lines

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,frames=len(x), interval=1, blit=True)
plt.show()
# anim.save('Test.mp4', writer=writer, dpi=400)