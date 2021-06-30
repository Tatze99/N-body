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

%matplotlib inline
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

def angle(x1,y1,z1,x2,y2,z2):
    return np.arccos((x1*x2+y1*y2+z1*z2)/(np.sqrt(x1**2+y1**2+z1**2)*np.sqrt(x2**2+y2**2+z2**2)))
    
def kinetic_energy(vx, vy, vz, m, n):
    energy = np.zeros(steps)
    for i in range(n):
        energy[:] += 0.5*m[i]*(vx[:,i]**2+vy[:,i]**2+vz[:,i]**2)
    return energy

def distance(x1,x2,y1,y2,z1,z2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def potential_energy(x, y, z, m, n):
    energy = np.zeros(steps)
    Matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            energy[:] += G*m[j]*m[i]/np.sqrt((x[:,i]-x[:,j])**2+(y[:,i]-y[:,j])**2+(z[:,i]-z[:,j])**2)
    return energy

def angular_momentum(x, y, z, vx, vy, vz, m, n):
    ang = np.zeros((steps, 3))
    for i in range(n):
        ang[:,0] += m[i]*(y[:,i]*vz[:,i]-z[:,i]*vy[:,i])
        ang[:,1] += m[i]*(z[:,i]*vx[:,i]-x[:,i]*vz[:,i])
        ang[:,2] += m[i]*(x[:,i]*vy[:,i]-y[:,i]*vx[:,i])
    return ang

def Laplace_Integral(x,y,z,vx,vy,vz,m,n):
    Laplace = np.zeros((steps, 3))
    cx = y*vz-z*vy
    cy = z*vx-x*vz
    cz = x*vy-y*vx
    r = np.sqrt(x**2+y**2+z**2)

    for i in range(n):
        Laplace[:,0] += cy[:,i]*vz[:,i]-cz[:,i]*vy[:,i]+G*m[i]*x[:,i]/r[:,i]
        Laplace[:,1] += cz[:,i]*vx[:,i]-cx[:,i]*vz[:,i]+G*m[i]*y[:,i]/r[:,i]
        Laplace[:,2] += cx[:,i]*vy[:,i]-cy[:,i]*vx[:,i]+G*m[i]*z[:,i]/r[:,i]
    return Laplace

#%%
# Initialize the data
# %matplotlib inline
command = "calc_t"
Input = np.loadtxt("Input.csv",delimiter=';') # input vx, vy, vz now in a.u. per year!!!!
Daten = np.loadtxt(command+"-solution.csv",delimiter=';')
mass = np.loadtxt("Input.csv",delimiter=';',usecols=[6])
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto', 'Sonde']

steps = len(Daten[:,0])
time  = Daten[:,0]

n = int((len(Daten[0,:])-1)/6)    # total number of planets
number = 11      # number of planets to display
if number > n: print("Error, too many planets to display")

# create the variables and assign them their values via a loop
var_names = ["x", "y", "z","vx", "vy", "vz"]
for i,name in enumerate(var_names):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]


# Plot the trajectories
%matplotlib auto
plt.figure(dpi=300)

# plt.figure(dpi=300, figsize=(2.5,3))
plt.plot(x[:,0:number], y[:,0:number],'.',markersize=1, label=Namen[0:number])
# plt.xlim(-33,70)
# plt.ylim(-40,40)
plt.xlim(-12,12)
plt.ylim(-12,12)
plt.legend(title='Planets')
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
# plt.yticks([])
# plt.title('Trajectories of the planets for 248 years')
# plt.savefig("Trajectories2D_"+command+"_Ausschnitt.pdf")
#%%
for i in range(n):
    print(np.sqrt(Input[i,0]**2+Input[i,1]**2+Input[i,2]**2))

#%%
# Calculate starting velocity of satellite
var_names2 = ["X", "Y", "Z","VX", "VY", "VZ"]
for i,name in enumerate(var_names2):
  globals()[name] = Input_tend[3,i]

ve = np.sqrt(VX**2+VY**2+VZ**2)
re = 4.2644e-5

xsat = X+VX*re/ve
ysat = Y+VY*re/ve
zsat = Z+VZ*re/ve

# vsat = 0
# for i in range(9): # without satellite
#     vsat += abs(mass[i]*(1/distance(xsat,x[0,i], ysat, y[0,i], zsat, z[0,i])-1/distance(x[0,i],x[0,9], y[0,i], y[0,9], z[0,i], z[0,9])))

# vsat += mass[9]/(2*r0)
# vsat = np.sqrt(vsat)*np.sqrt(2*G)

vsat = 9.170599

vxsat = VX*vsat/ve
vysat = VY*vsat/ve
vzsat = VZ*vsat/ve

probe_params = [xsat, ysat, zsat, vxsat, vysat, vzsat, 0, 0]
Input_tend[10,:] = probe_params
np.savetxt("Input_tend.csv", Input_tend, fmt='%1.20f', delimiter=';')
print(xsat, ysat, zsat, vxsat, vysat, vzsat)

#%%
Input_tend = np.loadtxt("Input_tend.csv",delimiter=';')
np.savetxt("Input_tend.csv", Input_tend, fmt='%1.20f', delimiter=';')
print(angle(Input_tend[3,0],Input_tend[3,1],Input_tend[3,2],Input_tend[5,0],Input_tend[5,1],Input_tend[5,2]))
print(Input_tend[3,0],Input_tend[3,1],Input_tend[3,2],Input_tend[5,0],Input_tend[5,1],Input_tend[5,2])

#%%
k = len(x[:,0])-1
j = 10
# Input_tend = Input
for i in range(j):
    Input_tend[i,0:6] = [x[k,i], y[k,i], z[k,i], vx[k,i], vy[k,i], vz[k,i]]
    
print(angle(x[k,3],y[k,3],z[k,3],x[k,5],y[k,5],z[k,5]))
# np.savetxt("Input_tend.csv", Input_tend, fmt='%1.20f', delimiter=';')

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
# Plot the angular momentum 
ang = angular_momentum(x,y,z,vx,vy,vz,mass,n)
Lap = Laplace_Integral(x,y,z,vx,vy,vz,mass,n)

plt.figure(dpi=300)
legend = ['$L_x$', '$L_y$', '$L_z$']
for i,energy in enumerate(legend):
    plt.plot(time,ang[:,i],label=legend[i])
plt.legend()

plt.figure(dpi=300)
legend = ['Laplacian $x$', 'Laplacian $y$', 'Laplacian $z$']
for i,energy in enumerate(legend):
    plt.plot(time,Lap[:,i],label=legend[i])
plt.legend()
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
=======
legend = ['Laplacian $x$', 'Laplacian $y$', 'Laplacian $z$']
for i,energy in enumerate(legend):
    plt.plot(time,Lap[:,i],label=legend[i])
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
%matplotlib auto
col = ['#000000', '#0C5DA5', '#0C5DA5', '#0C5DA5', '#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']
from matplotlib import animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=120, metadata=dict(artist='Martin Beyer'), bitrate=-1)

fig = plt.figure(figsize=(6,3))
ax1 = plt.axes()
line, = ax1.plot([], [], lw=2)
plt.xlabel('$x$ in a.u.')
plt.ylabel('$y$ in a.u.')

# for i in range(number):
#     ax1.plot(x[:,i],y[:,i],'-')  

lines = []
x2 = x[:,:number]
y2 = y[:,:number]
for i in range(number):
    lobj = ax1.plot([],[],'o', label=Namen[i], c=col[i])[0]
    lines.append(lobj)
plt.legend(loc='right')
plt.xlim(-6,9)
plt.ylim(-6,6)

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
anim = animation.FuncAnimation(fig, animate, init_func=init,frames=len(x), interval=1, blit=True)

plt.show()
# anim.save('Test.mp4', writer=writer, dpi=400)