"""
Created on Tue Apr 20th, 2021

Auswertung des N-body Problems
@author: Martin, Florian, Sebastian
"""
#%% header----------------------------------------
#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import rcParams
#rcParams['xtick.direction'] = 'in'
#rcParams['ytick.direction'] = 'in'
#rcParams['legend.fontsize'] = 7
#rcParams['font.size'] = 8
#rcParams['savefig.bbox'] = 'tight'
#rcParams['figure.figsize'] = (6,3)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.style.use(['science'])
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['legend.fontsize'] = 8
rcParams['font.size'] = 8
rcParams['savefig.bbox'] = 'tight'
rcParams['figure.figsize'] = (6,3)
rcParams['legend.frameon'] ='true'
rcParams['legend.framealpha'] = 0.74
plt.rcParams["font.family"] = "Times New Roman"
col = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e', '#e377c2', '#8c564b', '#17becf', '#bcbd22']

#%% Define functions----------------------------------------
#%matplotlib inline
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

def e_kin_nojup(vx, vy, vz, m, n):
    energy = np.zeros(steps)
    for i in range(n):
        if i == 5:
            continue
        else:
            energy[:] += 0.5*m[i]*(vx[:,i]**2+vy[:,i]**2+vz[:,i]**2)
    return energy

def e_pot_nojup(x, y, z, m, n):
    energy = np.zeros(steps)
    Matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            if i == 5:
                continue
            elif j == 5:
                continue
            else:
                energy[:] -= G*m[j]*m[i]/np.sqrt((x[:,i]-x[:,j])**2+(y[:,i]-y[:,j])**2+(z[:,i]-z[:,j])**2)
    return energy

def distance(x1,x2,y1,y2,z1,z2):
    return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def potential_energy(x, y, z, m, n):
    energy = np.zeros(steps)
    Matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            energy[:] -= G*m[j]*m[i]/np.sqrt((x[:,i]-x[:,j])**2+(y[:,i]-y[:,j])**2+(z[:,i]-z[:,j])**2)
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
%matplotlib inline
#command = "Swing"
command = "rk4"
Input = np.loadtxt("Input.csv",delimiter=';') # input vx, vy, vz now in a.u. per year!!!!
Daten = np.loadtxt(command+"-solution.csv",delimiter=';')
mass = np.loadtxt("Input.csv",delimiter=';',usecols=[6])
#Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto', 'Sonde']
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptun', 'Pluto']


#%%
# Plot in phase space
fig, ax= plt.subplots(3,dpi=400,figsize=(6,6))
ax[0].plot(x[:,0:number], vx[:,0:number],'-',lw=1, label='phase space x')
ax[1].plot(y[:,0:number], vy[:,0:number],'-',lw=1, label='phase space y')
ax[2].plot(z[:,0:number], vz[:,0:number],'-',lw=1, label={'Sun', 'Mercury', 'Venus'})
xlabel=['$x$  /  AU','$y$  /  AU','$z$  /  AU']
ylabel=['$v_x$  /  AU/yr','$v_y$  /  AU/yr','$v_z$  /  AU/yr']
for i in range(3):
    ax[i].set_xlabel(xlabel[i])
    ax[i].set_ylabel(ylabel[i])
    if (i!=2):
        ax[i].set_xticks([])
    ax[i].set_xlim(-32,45)
    ax[i].set_ylim(-13,10)
plt.legend(title='Planets', bbox_to_anchor=(1.26, 2.5))
#ax[0].set_title('Phase space plots')

#%%
# Plot the velocities against each other
fig, ax= plt.subplots(3,dpi=400,figsize=(6,6))
ax[0].plot(vx[:,0:number], vy[:,0:number],'-',lw=0.5, label=Namen[0:number])
ax[1].plot(vx[:,0:number], vz[:,0:number],'-',lw=0.5, label=Namen[0:number])
ax[2].plot(vy[:,0:number], vz[:,0:number],'-',lw=0.5, label=Namen[0:number])
xlabel=['$v_x$  /  AU/yr','$v_x$  /  AU/yr','$v_y$  /  AU/yr']
ylabel=['$v_y$  /  AU/yr','$v_z$  /  AU/yr','$v_z$  /  AU/yr']
for i in range(3):
    ax[i].set_xlabel(xlabel[i])
    ax[i].set_ylabel(ylabel[i])
    if (i!=2):
        ax[i].set_xticks([])
    ax[i].set_xlim(-13,11)
    ax[i].set_ylim(-10,10)
#plt.legend(title='Planets', bbox_to_anchor=(1.26, 2.5))
#ax[0].set_title('Phase space plots')

# plt.savefig("Phase_space_"+command+".pdf", bbox_inches='tight')

    