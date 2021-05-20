"""
Created on Tue Apr 20th, 2021

Auswertung des N-body Problems
@author: Martin, Florian, Sebastian
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 8
rcParams['savefig.bbox'] = 'tight'
rcParams['figure.figsize'] = (6,3)


G = 4*np.pi**2

def Schwerpunkt(x, y, z, m, n):
    xs = np.zeros(len(x[:,0]))
    ys = np.zeros(len(x[:,0]))
    zs = np.zeros(len(x[:,0]))
    
    for i in range(n):
        xs[:] += m[i]*x[:,i]
        ys[:] += m[i]*x[:,i]
        zs[:] += m[i]*x[:,i]
    return xs, ys, zs

def kinetic_energy(vx, vy, vz, m, n):
    energy = np.zeros(steps)
    for i in range(n):
        energy[:] += 0.5*m[i]*(vx[:,i]**2+vy[:,i]**2+vz[:,i]**2)
    return energy

def distance(x1,x2,y1,y2,z1,z2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def potential_energy(x, y, z, m, n):
    energy = np.zeros(steps)
    for i in range(n):
        for j in range(i):
            energy[:] -= G*m[j]*m[i]/np.sqrt((x[:,i]-x[:,j])**2+(y[:,i]-y[:,j])**2+(z[:,i]-z[:,j])**2)
    return energy

command = "rk4"
Daten = np.loadtxt(command+"-solution.csv",delimiter=';')
mass = np.loadtxt("Input.csv",delimiter=';',usecols=[6])
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto']

steps = len(Daten[:,0])
time  = Daten[:,0]
number = 10        # number of planets to display
n = len(mass)-1     # total number of planets

# create the variables and assign them their values via a loop
var_names = ["x", "y", "z","vx", "vy", "vz"]
for i,name in enumerate(var_names):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]

#%%
# Calculate starting velocity of satellite
r0 = 6.685e-6
ve = np.sqrt(vx[0,3]**2+vy[0,3]**2+vz[0,3]**2)
re = 4.2644e-5

xsat = x[0,3]+vx[0,3]*re/ve
ysat = y[0,3]+vy[0,3]*re/ve
zsat = z[0,3]+vz[0,3]*re/ve

vsat = 0
for i in range(n-1):
    vsat += abs(mass[i]*(1/distance(xsat,x[0,i], ysat, y[0,i], zsat, z[0,i])-1/distance(x[0,i],x[0,9], y[0,i], y[0,9], z[0,i], z[0,9])))

vsat += mass[9]/(2*r0)
vsat = np.sqrt(vsat)*np.sqrt(2*G)
print(vsat)
#%%
# Plot the trajectories
plt.figure(dpi=300)
# plt.figure(dpi=300, figsize=(2.5,3))
plt.plot(x[:,0:number], y[:,0:number],'.',markersize=1, label=Namen[0:number])
plt.xlim(-33,70)
# plt.ylim(-3.99,3.99)
# plt.xlim(-4,4)
plt.legend(title='Planets')
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
# plt.yticks([])
# plt.title('Trajectories of the planets for 248 years')
# plt.savefig("Trajectories2D_"+command+"_Ausschnitt.pdf")

#%%
# Plot in phase space
fig, ax= plt.subplots(3,dpi=200,figsize=(6,6))
ax[0].plot(x[:,0:number], vx[:,0:number],'-',lw=1, label='phase space x')
ax[1].plot(y[:,0:number], vy[:,0:number],'-',lw=1, label=Namen[0:number])
ax[2].plot(z[:,0:number], vz[:,0:number],'-',lw=1, label=Namen[0:number])
xlabel=['x','y','z']
ylabel=['vx','vy','vz']
for i in range(3):
    ax[i].set_xlabel(xlabel[i])
    ax[i].set_ylabel(ylabel[i])
    if (i!=2):
        ax[i].set_xticks([])
    ax[i].set_xlim(-40,40)
    ax[i].set_ylim(-8,8)
plt.legend(title='Planets', bbox_to_anchor=(1.26, 2.5))
ax[0].set_title('Phase space plots')

# plt.savefig("Phase_space_"+command+".pdf", bbox_inches='tight')
#%%
# Plot the energy
energy_conversion = 4.47 # e37
kin = energy_conversion*kinetic_energy(vx, vy, vz, mass, number)
pot = energy_conversion*potential_energy(x, y, z, mass, number)

plt.figure(dpi=300)
legend = ['kinetic energy', 'potential energy', 'total energy']
for i,energy in enumerate([kin, pot, kin+pot]):
    plt.plot(time,energy,label=legend[i])
plt.legend()
plt.xlabel('time $t$ in years')
plt.ylabel('energy $E$ in $10^37$ J')
# plt.title('Energy of the system as a function of time')
plt.savefig("Energy_"+command+".pdf")

# print the kinetic energies of the planets
# for i in range(n):
#     print(0.5*mass[i]*(vx[-1,i]**2+vy[-1,i]**2+vz[-1,i]**2))
#%%
# Plot in 3D
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(projection='3d')

for i,name in enumerate(Namen):
    ax.plot(x[:,i], y[:,i], z[:,i], label=Namen[i])

# ax.set_title('Trajectories of all planets')
ax.set_xlabel('$x$ in AU')
ax.set_ylabel('$y$ in AU')
ax.set_zlabel('$z$ in AU')
ax.legend(loc='center left')
# plt.savefig("Trajectories3D_"+command+".pdf")

#%%
# Plot Schwerpunkt
xs, ys, zs = Schwerpunkt(x,y,z,mass, number)
print(xs[0], ys[0], zs[0])
print(xs[-1], ys[-1], zs[-1])
#%%
from matplotlib import animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=120, metadata=dict(artist='Martin Beyer'), bitrate=-1)

fig = plt.figure(figsize=(6,3))
ax1 = plt.axes()
line, = ax1.plot([], [], lw=2)
plt.xlabel('$x$ in a.u.')
plt.ylabel('$y$ in a.u.')

for i in range(number):
    ax1.plot(x[:,i],y[:,i],'-')  

lines = []
x2 = x[:,:number]
y2 = y[:,:number]
for index in range(number):
    lobj = ax1.plot(x[:1,index],y[:1,index],'o', label=Namen[index])[0]
    lines.append(lobj)
plt.legend(loc='right')
# plt.xlim(-35,70)

def init():
    for line in lines:
        line.set_data([],[])
    return lines

def animate(i):
    for lnum,line in enumerate(lines):
        # index=-i # Standard
        index = -1
        # index = round(-0.012*lnum**6)
        if(i+index > 0):
            line.set_data(x[i+index:i,lnum], y[i+index:i,lnum]) # set data for each line separately.
        else:
            line.set_data(x[:i,lnum], y[:i,lnum])

    return lines

# call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,frames=len(x), interval=1, blit=True)
plt.show()
# anim.save('Test.mp4', writer=writer, dpi=400)