"""
Created on Tue Apr 20th, 2021

Auswertung des N-body Problems
@author: Martin, Florian, Sebastian
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

plt.style.use(['science'])  ## pip install SciencePlots
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['legend.fontsize'] = 8
rcParams['font.size'] = 8
rcParams['savefig.bbox'] = 'tight'
rcParams['figure.figsize'] = (6,3)
rcParams['legend.frameon'] ='true'
rcParams['legend.framealpha'] = 0.74
plt.rcParams["font.family"] = "Times New Roman"


PDF = Path("pdf_files")
CSV = Path("csv_files")

%matplotlib inline
G = 4*np.pi**2
col = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e', '#e377c2', '#8c564b', '#17becf', '#bcbd22']

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

def potential_energy(x, y, z, m, n):
    energy = np.zeros(steps)
    Matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            energy[:] += G*m[j]*m[i]/np.sqrt((x[:,i]-x[:,j])**2+(y[:,i]-y[:,j])**2+(z[:,i]-z[:,j])**2)
    return energy

def total_energy(Daten,m,n):
    x = Daten[:,0*n+1:(0+1)*n+1]
    y = Daten[:,1*n+1:(1+1)*n+1]  
    z = Daten[:,2*n+1:(2+1)*n+1]  
    vx = Daten[:,3*n+1:(3+1)*n+1]  
    vy = Daten[:,4*n+1:(4+1)*n+1]  
    vz = Daten[:,5*n+1:(5+1)*n+1]  
    return kinetic_energy(vx, vy, vz, m, n)-potential_energy(x, y, z, m, n)
    
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

def Data_to_variables(Daten, n):
    t = Daten[:,0]
    x = Daten[:,0*n+1:(0+1)*n+1]
    y = Daten[:,1*n+1:(1+1)*n+1]  
    z = Daten[:,2*n+1:(2+1)*n+1]  
    vx = Daten[:,3*n+1:(3+1)*n+1]  
    vy = Daten[:,4*n+1:(4+1)*n+1]  
    vz = Daten[:,5*n+1:(5+1)*n+1]
    return t,x,y,z,vx,vy,vz
#%%  Initialize the data and make a first plot
%matplotlib inline
Input = np.loadtxt(CSV/"Input.csv",delimiter=';') # input vx, vy, vz now in a.u. per year!!!!
Daten = np.loadtxt(CSV/"rk4-solution_248.csv",delimiter=';')
mass = np.loadtxt(CSV/"Input.csv",delimiter=';',usecols=[6])
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'probe', 'probe2', 'probe3']

steps = len(Daten[:,0])
time  = Daten[:,0]

n = int((len(Daten[0,:])-1)/6)    # total number of planets
number = n                        # number of planets to display

# create the variables and assign them their values via a loop
for i,name in enumerate(["x", "y", "z","vx", "vy", "vz"]):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]


# Plot the trajectories
plt.figure(dpi=300)
for i in range(number):
    plt.plot(x[:,i], y[:,i], label=Namen[i], c=col[number-1-i])
plt.xlim(-33,70)
plt.ylim(-40,40)
# plt.xlim(-5,5)
# plt.ylim(-5,5)

plt.legend()
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')

# plt.savefig(PDF/"Trajectories2D_"+command+"_Ausschnitt.pdf")


#%% Compare trajectories for rk4, lf, fwd

Daten_rk4 = np.loadtxt(CSV/"rk4-solution_248.csv",delimiter=';')
Daten_fwd = np.loadtxt(CSV/"fwd-solution_248.csv",delimiter=';')
Daten_lf = np.loadtxt(CSV/"lf-solution_248.csv",delimiter=';')
plt.plot(Daten_rk4[::9,0*n+2], Daten_rk4[::9,1*n+2],'.',markersize=1, label="Runge-Kutta", c=col[0])
plt.plot(Daten_fwd[:,0*10+2], Daten_fwd[:,1*10+2],'.',markersize=0.5, label="Euler-Forward", c=col[1])
plt.plot(Daten_lf[::9,0*n+2], Daten_lf[::9,1*n+2],'.',markersize=0.5, label="Leap-Frog", alpha=0.7, c=col[2])
plt.xlim(-1,1)
plt.ylim(-1,1)

plt.legend(title='Integrators', markerscale=5.)
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
plt.savefig(PDF/"Comparison_different_integrators.pdf")

#%%  Compare different error-values for Cash-Karp

t_18, x_18, y_18, z_18, vx_18, vy_18, vz_18 = Data_to_variables(np.loadtxt(CSV/"rk4-solution_e18.csv",delimiter=';'), n)
t_10, x_10, y_10, z_10, vx_10, vy_10, vz_10 = Data_to_variables(np.loadtxt(CSV/"rk4-solution_e10.csv",delimiter=';'), n)
t_7, x_7, y_7, z_7, vx_7, vy_7, vz_7 = Data_to_variables(np.loadtxt(CSV/"rk4-solution_e7.csv",delimiter=';'), n)
plt.plot(x_7[:,1],y_7[:,1],'.',markersize=1, label="$\\Delta_0 = 10^{-7}$")
plt.plot(x_10[:,1],y_10[:,1],'.',markersize=1, label="$\\Delta_0 = 10^{-10}$")
plt.plot(x_18[:,1],y_18[:,1],'.',markersize=1, label="$\\Delta_0 = 10^{-18}$")
plt.legend(title='Error specification', markerscale=5.)
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
# plt.savefig(PDF/"Comparison_different_errors.pdf")

#%% All trajectories calculated by angle + calc_sat for Jupiter, Saturn, Uranus, Neptun, Pluto
t_J,x_J, y_J, z_J, vx_J, vy_J, vz_J = Data_to_variables(np.loadtxt(CSV/"calc_t-solution_Jupiter.csv",delimiter=';'), n)
t_S, x_S, y_S, z_S, vx_S, vy_S, vz_S = Data_to_variables(np.loadtxt(CSV/"calc_t-solution_Saturn.csv",delimiter=';'), n)
t_U, x_U, y_U, z_U, vx_U, vy_U, vz_U = Data_to_variables(np.loadtxt(CSV/"calc_t-solution_Uranus.csv",delimiter=';'), n)
t_N, x_N, y_N, z_N, vx_N, vy_N, vz_N = Data_to_variables(np.loadtxt(CSV/"calc_t-solution_Neptun.csv",delimiter=';'), n)
t_P, x_P, y_P, z_P, vx_P, vy_P, vz_P = Data_to_variables(np.loadtxt(CSV/"calc_t-solution_Pluto.csv",delimiter=';'), n)

np.savetxt(CSV/"calc_t_Jupiter.csv", np.stack((x_J[::10,5],x_J[::10,10], y_J[::10,5],y_J[::10,10])), fmt='%.15f')
np.savetxt(CSV/"calc_t_Saturn.csv", np.stack((x_S[:,6],x_S[:,10], y_S[:,6],y_S[:,10])), fmt='%.15f')
np.savetxt(CSV/"calc_t_Uranus.csv", np.stack((x_U[::50,7],x_U[::50,10], y_U[::50,7],y_U[::50,10])), fmt='%.15f')
np.savetxt(CSV/"calc_t_Neptun.csv", np.stack((x_N[::100,8],x_N[::100,10], y_N[::100,8],y_N[::100,10])), fmt='%.15f')
np.savetxt(CSV/"calc_t_Pluto.csv", np.stack((x_P[::150,9],x_P[::150,10], y_P[::150,9],y_P[::150,10])), fmt='%.15f')

#%% Darstellung der Trajektorien für Sonden, die zu den Planeten fliegen nach der Berechnung von calc_angle + calc_t
Daten = np.loadtxt(CSV/"rk4-solution_248.csv",delimiter=';')
n = 10
for i,name in enumerate(["x", "y", "z","vx", "vy", "vz"]):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]
  
Jupiter = np.loadtxt(CSV/"calc_t_Jupiter.csv")
Saturn = np.loadtxt(CSV/"calc_t_Saturn.csv")
Uranus = np.loadtxt(CSV/"calc_t_Uranus.csv")
Neptun = np.loadtxt(CSV/"calc_t_Neptun.csv")
Pluto = np.loadtxt(CSV/"calc_t_Pluto.csv")

plt.figure()
Namen = ["Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

for i in range(4,9):
    plt.plot(x[::10,i], y[::10,i], c=col[9-i], alpha=0.2)

plt.plot(x[::100,3], y[::100,3], c=col[6])

plt.plot(Jupiter[0,:], Jupiter[2,:], label=Namen[0], c=col[4])
plt.plot(Jupiter[1,:], Jupiter[3,:], '--', c=col[4])

plt.plot(Saturn[0,:], Saturn[2,:], label=Namen[1], c=col[3])
plt.plot(Saturn[1,:], Saturn[3,:], '--', c=col[3])

plt.plot(Uranus[0,:], Uranus[2,:], label=Namen[2], c=col[2])
plt.plot(Uranus[1,:], Uranus[3,:], '--', c=col[2])

plt.plot(Neptun[0,:], Neptun[2,:], label=Namen[3], c=col[1])
plt.plot(Neptun[1,:], Neptun[3,:], '--', c=col[1])

plt.plot(Pluto[0,:], Pluto[2,:], label=Namen[4], c=col[0])
plt.plot(Pluto[1,:], Pluto[3,:], '--', c=col[0])

plt.legend(title='Planets')
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
plt.xlim(-33,20)
plt.ylim(-33,10)
# plt.savefig(PDF/"Trajectories_calc_t.pdf")


#%% Trajectories of the New Horizon space probe (original trajectory and adjusted trajectories by the algorithm "Swing")

Daten = np.loadtxt(CSV/"NewHorizon-solution.csv",delimiter=';')
n = int((len(Daten[0,:])-1)/6)
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
col = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e', '#e377c2', '#8c564b', '#17becf', '#bcbd22']
# col = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e', '#e377c2', '#8c564b', '#17becf', '#bcbd22']
for i,name in enumerate(["x", "y", "z","vx", "vy", "vz"]):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]
time  = Daten[:,0]

# Plot the trajectories
plt.figure(dpi=300)

# plt.figure(dpi=300, figsize=(2.5,3))
for index,i in enumerate([3,5,9]):
    plt.plot(x[:,i], y[:,i], label=Namen[i], c=col[2-index])

for index,i in enumerate([3,5,9]):
    plt.plot(x[-1,i], y[-1,i],'.', c=col[2-index])
plt.plot(x[:,10], y[:,10], label="New Horizons (original data)", c=col[4])
plt.plot(x[:,11], y[:,11], '--', label="New Horizons (adjusted)", c=col[4])
plt.plot(x[:,12], y[:,12], ':', label="New Horizons (adjusted)", c=col[4])

# uncomment for zoom out
# plt.ylim(-35,8)
# plt.legend()
# plt.xlabel('$x$ in AU')
# plt.ylabel('$y$ in AU')
# plt.savefig(PDF/"NewHorizons_trajectories.pdf")

# uncomment for zoom in
plt.plot(x[2690,10], y[2690,10], 'o', c=col[4])
plt.plot(x[2690,11], y[2690,11], 'o', c=col[4])
plt.plot(x[2690,12], y[2690,12], 'o', c=col[4])
plt.plot(x[2690,5], y[2690,5], 'o', c=col[1])
plt.xlim(-1.95,-1.75)
plt.ylim(-4.8,-4.4)
plt.yticks([-4.5,-4.6,-4.7])
plt.xticks([-1.9,-1.8])
# plt.savefig(PDF/"NewHorizons_trajectories_zoom.pdf")


#%%  Plot the velocity of New Horizons as a function of time (zoom in and zoom out)
vel1 = np.sqrt(vx[:,10]**2+vy[:,10]**2+vz[:,10]**2)
vel2 = np.sqrt(vx[:,11]**2+vy[:,11]**2+vz[:,11]**2)
vel3 = np.sqrt(vx[:,12]**2+vy[:,12]**2+vz[:,12]**2)
plt.figure()
plt.plot(time, vel1, c=col[4], label="New Horizons (original data)")
plt.plot(time,vel2, '--', c=col[4])
plt.plot(time,vel3,':', c=col[4])

plt.xlabel('$t$ in yr')
plt.ylabel('velocity $v$ in AU per yr')
plt.legend()
# plt.savefig(PDF/"NewHorizons_Energies.pdf")

# plt.xlim(1.1,1.3)
# plt.ylim(3,5)
# plt.yticks([3.5,4,4.5])
# plt.xticks([1.15,1.2,1.25])
# plt.savefig(PDF/"NewHorizons_Energies_zoom.pdf")

#%% All trajectories calculated by Swing-by for Jupiter, Saturn, Uranus, Neptun, Pluto
n=11
t_J, x_J, y_J, z_J, vx_J, vy_J, vz_J = Data_to_variables(np.loadtxt(CSV/"Swing-by-solution_Jupiter.csv",delimiter=';'), n)
t_S, x_S, y_S, z_S, vx_S, vy_S, vz_S = Data_to_variables(np.loadtxt(CSV/"Swing-by-solution_Saturn.csv",delimiter=';'), n)
t_U, x_U, y_U, z_U, vx_U, vy_U, vz_U = Data_to_variables(np.loadtxt(CSV/"Swing-by-solution_Uranus.csv",delimiter=';'), n)
t_N, x_N, y_N, z_N, vx_N, vy_N, vz_N = Data_to_variables(np.loadtxt(CSV/"Swing-by-solution_Neptun.csv",delimiter=';'), n)
t_P, x_P, y_P, z_P, vx_P, vy_P, vz_P = Data_to_variables(np.loadtxt(CSV/"calc_t-solution_Pluto.csv",delimiter=';'), n)

np.savetxt(CSV/"Swing-by_Jupiter.csv", np.stack((x_J[:,5],x_J[:,10], y_J[:,5],y_J[:,10])), fmt='%.15f')
np.savetxt(CSV/"Swing-by_Saturn.csv", np.stack((x_S[:,6],x_S[:,10], y_S[:,6],y_S[:,10])), fmt='%.15f')
np.savetxt(CSV/"Swing-by_Uranus.csv", np.stack((x_U[:,7],x_U[:,10], y_U[:,7],y_U[:,10])), fmt='%.15f')
np.savetxt(CSV/"Swing-by_Neptun.csv", np.stack((x_N[:,8],x_N[:,10], y_N[:,8],y_N[:,10])), fmt='%.15f')
np.savetxt(CSV/"Swing-by_Pluto.csv", np.stack((x_P[:,9],x_P[:,10], y_P[:,9],y_P[:,10])), fmt='%.15f')

#%% First evaulate previous cell!
Daten = np.loadtxt(CSV/"rk4-solution_248.csv",delimiter=';')
time  = Daten[:,0]
n = 10
for i,name in enumerate(["x", "y", "z","vx", "vy", "vz"]):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]
  
Jupiter = np.loadtxt(CSV/"Swing-by_Jupiter.csv")
Saturn = np.loadtxt(CSV/"Swing-by_Saturn.csv")
Uranus = np.loadtxt(CSV/"Swing-by_Uranus.csv")
Neptun = np.loadtxt(CSV/"Swing-by_Neptun.csv")
Pluto = np.loadtxt(CSV/"Swing-by_Pluto.csv")

plt.figure()
Namen = ["Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

for i in range(4,8):
    plt.plot(x[::10,i], y[::10,i], c=col[9-i], alpha=0.2)

plt.plot(x[::100,3], y[::100,3], c=col[6])

plt.plot(Jupiter[0,:], Jupiter[2,:], label=Namen[0], c=col[4])
plt.plot(Jupiter[1,:], Jupiter[3,:], '--', c=col[4])

plt.plot(Saturn[0,:], Saturn[2,:], label=Namen[1], c=col[3])
plt.plot(Saturn[1,:], Saturn[3,:], '--', c=col[3])

plt.plot(Uranus[0,:], Uranus[2,:], label=Namen[2], c=col[2])
plt.plot(Uranus[1,:], Uranus[3,:], '--', c=col[2])

plt.plot(Neptun[0,:], Neptun[2,:], label=Namen[3], c=col[1])
plt.plot(Neptun[1,:], Neptun[3,:], '--', c=col[1])

plt.plot(Pluto[0,:], Pluto[2,:], label=Namen[4], c=col[0])
plt.plot(Pluto[1,:], Pluto[3,:], '--', c=col[0])

plt.legend(title='Planets')
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
plt.xlim(-30,30)
plt.ylim(-35,20)

# plt.savefig(PDF/"Swing-by_trajectories.pdf")
#%%
Nam = ["Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
vel1 = np.sqrt(vx_J[:,10]**2+vy_J[:,10]**2+vz_J[:,10]**2)
vel2 = np.sqrt(vx_S[:,10]**2+vy_S[:,10]**2+vz_S[:,10]**2)
vel3 = np.sqrt(vx_U[:,10]**2+vy_U[:,10]**2+vz_U[:,10]**2)
vel4 = np.sqrt(vx_N[:,10]**2+vy_N[:,10]**2+vz_N[:,10]**2)
vel5 = np.sqrt(vx_P[:,10]**2+vy_P[:,10]**2+vz_P[:,10]**2)
plt.figure()
    
plt.axhline(y=1.4, color='gray', lw=0.7, alpha=0.3)
plt.axhline(y=4, color='gray', lw=0.7, alpha=0.3)
plt.plot(t_J, vel1, c=col[4], label="Jupiter")
plt.plot(t_S, vel2, c=col[3], label="Saturn")
plt.plot(t_U, vel3, c=col[2], label="Uranus")
plt.plot(t_N, vel4, c=col[1], label="Neptune")
plt.plot(t_P, vel5, c=col[0], label="Pluto")

plt.xlim(-0.5,10)

plt.xlabel('$t$ in yr')
plt.ylabel('velocity $v$ in AU per yr')
plt.legend()
plt.savefig(PDF/"Swing-by_Energies.pdf")
#%%
# Calculate starting velocity of satellite
var_Names2 = ["X", "Y", "Z","VX", "VY", "VZ"]
for i,name in enumerate(var_Names2):
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
np.savetxt(CSV/"Input_tend.csv", Input_tend, fmt='%1.20f', delimiter=';')
print(xsat, ysat, zsat, vxsat, vysat, vzsat)

#%%  Cut big files in half by deleting every second time step
Daten2 = np.loadtxt(CSV/"Swing-by_Pluto.csv",delimiter=';')
Daten2 = np.delete(Daten2, (np.arange(0,len(Daten2[:,0]),2)), axis=0)

n=13
for i,name in enumerate(["x", "y", "z","vx", "vy", "vz"]):
  globals()[name] = Daten2[:,i*n+1:(i+1)*n+1]

plt.figure(dpi=300)

plt.plot(x[:,1:number], y[:,1:number],'.',markersize=1)
# np.savetxt(CSV/"Swing-by_Pluto.csv", Daten2, fmt='%1.15f', delimiter=';')
#%%  Plot in phase space
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

# plt.savefig(PDF/"Phase_space_"+command+".pdf")
#%%  Plot the energy
energy_conversion = 4.47 # e37
kin = energy_conversion*kinetic_energy(vx, vy, vz, mass, number)
pot = energy_conversion*potential_energy(x, y, z, mass, number)

plt.figure(dpi=300)
legend = ['kinetic energy', 'potential energy', 'total energy']
for i,energy in enumerate([kin, pot, kin-pot]):
    plt.plot(time,energy,label=legend[i])
plt.legend()
plt.xlabel('time $t$ in years')
plt.ylabel('energy $E$ in $10^37$ J')
# plt.title('Energy of the system as a function of time')
plt.savefig(PDF/"Energy_"+command+".pdf")

#%% Compare total energy for rk4, lf, fwd
E_rk4 = 4.47*total_energy(Daten_rk4,mass,n)
E_fwd = 4.47*total_energy(Daten_fwd,mass,n)
E_lf = 4.47*total_energy(Daten_lf,mass,n)
plt.plot(time,E_rk4, label="Runge-Kutta")
plt.plot(time,E_fwd, label="Euler-forward")
plt.plot(time,E_lf, label="leap frog")
plt.legend(title="Integrators")
plt.xlabel('time $t$ in years')
plt.ylabel('energy $E$ in $10^{37}$ J')

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
#%% Plot in 3D
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(projection='3d')

for i in range(n):
    ax.plot(x[:,i], y[:,i], z[:,i], label=Namen[i], c=col[n-1-i])

# ax.set_title('Trajectories of all planets')
ax.set_xlabel('$x$ in AU')
ax.set_ylabel('$y$ in AU')
ax.set_zlabel('$z$ in AU')
ax.legend(loc='center left')
ax.minorticks_off()

# plt.savefig(PDF/"Trajectories3D_"+command+".pdf")

#%% Plot Schwerpunkt (center of gravity)
xs, ys, zs = Schwerpunkt(x,y,z,mass, number)
print(xs[0], ys[0], zs[0])
print(xs[-1], ys[-1], zs[-1])

legend = ['Laplacian $x$', 'Laplacian $y$', 'Laplacian $z$']
for i,energy in enumerate(legend):
    plt.plot(time,Lap[:,i],label=legend[i])
plt.legend()


#%% Plot probe velocity to minimum distance
v_dep_Mars = np.loadtxt(CSV/"velocity_dependence_Mars.csv",delimiter=';')
v_dep_Jupiter = np.loadtxt(CSV/"velocity_dependence_Jupiter.csv",delimiter=';')
v_dep_Saturn = np.loadtxt(CSV/"velocity_dependence_Saturn.csv",delimiter=';')
v_dep_Uranus = np.loadtxt(CSV/"velocity_dependence_Uranus.csv",delimiter=';')
v_dep_Neptune = np.loadtxt(CSV/"velocity_dependence_Neptune.csv",delimiter=';')
v_dep_Pluto = np.loadtxt(CSV/"velocity_dependence_Pluto.csv",delimiter=';')

plt.rcParams["figure.figsize"] = (6,3)
plt.figure()
plt.plot(v_dep_Mars[:,0], v_dep_Mars[:,1], label="Mars")
plt.plot(v_dep_Jupiter[:,0], v_dep_Jupiter[:,1], label="Jupiter")
plt.plot(v_dep_Saturn[:,0], v_dep_Saturn[:,1], label="Saturn")
plt.plot(v_dep_Uranus[:,0], v_dep_Uranus[:,1], label="Uranus")
plt.plot(v_dep_Neptune[:,0], v_dep_Neptune[:,1], label="Neptune")
plt.plot(v_dep_Pluto[:,0], v_dep_Pluto[:,1], label="Pluto")
plt.xlim(9.28, 9.92)
plt.ylim(-0.4,8)
plt.xlabel("velocity $v$ AU/yr")
plt.ylabel("minimum distance in AU")
plt.legend()
plt.annotate("0.00063", (9.4,-0.15), c=col[1])
plt.plot(v_dep_Jupiter[44,0], v_dep_Jupiter[44,1], 'o', c=col[1])
plt.annotate("0.00064", (9.35,1), c=col[2])
plt.plot(v_dep_Saturn[27,0], v_dep_Saturn[27,1], 'o', c=col[2])
plt.annotate("0.139", (9.63,0.139), c=col[3])
plt.plot(v_dep_Uranus[27,0], v_dep_Uranus[27,1], 'o', c=col[3])
plt.annotate("0.358", (9.76,0.358), c=col[4])
plt.plot(v_dep_Neptune[32,0], v_dep_Neptune[32,1], 'o', c=col[4])
plt.annotate("1.345", (9.78,1.345), c=col[5])
plt.plot(v_dep_Pluto[25,0], v_dep_Pluto[25,1], 'o', c=col[5])
plt.savefig(PDF/"velocity_dependence.pdf")

# plt.plot(v_dep_Mars[:,0], v_dep_Mars[:,1], label="Mars")
# plt.yticks([0.04, 0.08])
# plt.xticks([8.6, 8.7])
# plt.annotate("0.0014", (8.61,0.014), c=col[0])
# plt.plot(v_dep_Mars[42,0], v_dep_Mars[42,1], 'o', c=col[0])
# plt.savefig(PDF/"velocity_dependence_Zoom.pdf")
#%%
#Load data for convergence plots
#Daten = np.loadtxt(command+"-convergence.csv",delimiter=';')
Daten = np.loadtxt(CSV/"Convergence/rk4-cash-karp-convergence-v.csv",delimiter=';')
h,f = Daten[:,0], Daten[:,1]

#For Euler convergence
# stepsize = 2*np.logspace(-8,-1, num=16)
# Error = stepsize
#Fit = 0.4438*(stepsize**0.99536)
# Fit = 0.41736*(stepsize**0.99894)

#For RK4 convergence
# stepsize = 4*np.logspace(-4,-1, num=16)
# Error = stepsize**4
#Fit = 10.5E-4 * (stepsize**3.9640)
# Fit = 11.14E-5 * (stepsize**3.9132)

#For RK5 convergence
# stepsize = 4*np.logspace(-3,-1, num=16)
# Error = stepsize**5
#Fit = 14.94E-6 * (stepsize**4.8168)
# Fit = 5.779E-6 * (stepsize**4.9308)

#For convergence of RK4-part of Cash-Karp
stepsize = 7*np.logspace(-4,-1, num=16)
Error = stepsize**4
# Fit = 111.E-6 * (stepsize**3.9712)
Fit = 8.13E-6 * (stepsize**3.8864)

#For lf convergence
# stepsize = 4*np.logspace(-7,-1, num=16)
# Error = stepsize**2
# Fit = 16.7E-3 * (stepsize**1.9773)
#Fit = 4.708E-6 * (stepsize**4.832)

# Convergence plot
rcParams['figure.figsize'] = (6,3)
plt.figure(dpi=400)
plt.loglog(h, f, 'o-', label='Behaviour of integrator')
plt.loglog(stepsize, Error, '--',label='$f(x)=x^4$')

plt.loglog(stepsize, Fit, '--',label='Fit: $f(x)=8.13\\cdot 10^{-6}\\cdot x^{3.8864}$') # CK RK5
# plt.loglog(stepsize, Fit, '--',label='Fit: $f(x)=5.78\\cdot 10^{-6}\\cdot x^{4.9308}$') # CK RK5
# plt.loglog(stepsize, Fit, '--',label='Fit: $f(x)=1.67\\cdot 10^{-2}\\cdot x^{1.9773}$') # leap frog
# plt.loglog(stepsize, Fit, '--',label='Fit: $f(x)=0.417\\cdot x^{0.9989}$') # euler velocity
# plt.loglog(stepsize, Fit, '--',label='Fit: $f(x)=111\\cdot 10^{-6}\\cdot x^{3.9132}$') # RK4 velocity
plt.legend(prop={'size': 9})
plt.xlabel('stepsize $h$')
# plt.ylabel('Error in positions')
plt.ylabel('Error in velocities')
# plt.savefig(PDF/"RK4-Cash-Karp-convergence-v.pdf")
# plt.savefig(PDF/"rk5 convergence-v.pdf")
# plt.savefig(PDF/"lf convergence-v.pdf")
# plt.savefig(PDF/"rk4 convergence-v.pdf")

#%%  Plot the solar system with different satellite trajectories (vmin, vmax for every planet)
Daten = np.loadtxt(CSV/"sat-trajectories-solution.csv",delimiter=';')
Namen = ['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Mars-Sonde','Mars-Sonde','Jupiter-Sonde','Saturn-Sonde','Uranus-Sonde','Neptun-Sonde','Pluto-Sonde']
mass = np.loadtxt(CSV/"Input.csv",delimiter=';',usecols=[6])

steps = len(Daten[:,0])
time  = Daten[:,0]

n = int((len(Daten[0,:])-1)/6)    # total number of planets
number = 10
#number = 11
if number > n: print("Error, too many planets to display")

# create the variables and assign them their values via a loop
var_names = ["x", "y", "z","vx", "vy", "vz"]
for i,name in enumerate(var_names):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]

# Plot the trajectories
plt.figure(dpi=400)
for i in range(number):
    plt.plot(x[:,i], y[:,i],'-',markersize=0.5, label=Namen[i], c=col[number-1-i])

plt.plot(x[:130,11], y[:130,11],'--',markersize=0.5, c=col[4])
plt.plot(x[:285,12], y[:285,12],'--',markersize=0.5, c=col[3])
plt.plot(x[:750,13], y[:750,13],'--',markersize=0.5, c=col[2])
plt.plot(x[:1500,14], y[:1500,14],'--',markersize=0.5, c=col[1])
plt.plot(x[:2660,15], y[:2660,15],'--',markersize=0.5, c=col[0])

Daten = np.loadtxt(CSV/"sat-trajectories-solution-min.csv",delimiter=';')
var_names = ["x", "y", "z","vx", "vy", "vz"]
for i,name in enumerate(var_names):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]
  
plt.plot(x[:130,11], y[:130,11],':',markersize=0.5, c=col[4])
plt.plot(x[:285,12], y[:285,12],':',markersize=0.5, c=col[3])
plt.plot(x[:750,13], y[:750,13],':',markersize=0.5, c=col[2])
plt.plot(x[:1500,14], y[:1500,14],':',markersize=0.5, c=col[1])
plt.plot(x[:1300,15], y[:1300,15],':',markersize=0.5, c=col[0])


plt.xlim(-50,45)
plt.ylim(-20,50)
#plt.xlim(-6,6)
#plt.ylim(-5,5)
plt.legend()
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
# plt.savefig(PDF/"Vmin_vmax_Planets.pdf")

#%% Energy of the solar system as a function of time with and without Jupiter

Daten = np.loadtxt(CSV/"rk4-solution_248.csv",delimiter=';')
for i,name in enumerate(["x", "y", "z","vx", "vy", "vz"]):
  globals()[name] = Daten[:,i*n+1:(i+1)*n+1]

# Plot the energy
energy_conversion = 4.47 # e37
#Florian: multiplied by 100 and set 10^35 for nicer y-labelling
kin = energy_conversion*kinetic_energy(vx, vy, vz, mass, number)*100
pot = energy_conversion*potential_energy(x, y, z, mass, number)*100

kin2 = energy_conversion*e_kin_nojup(vx, vy, vz, mass, number)*1000
pot2 = energy_conversion*e_pot_nojup(x, y, z, mass, number)*1000

plt.figure(dpi=400)
legende = ['$E_\\text{kin}$', '$E_\\text{pot}$', '$E_\\text{total}$']
for i,energy in enumerate([kin, pot, kin+pot]):
    plt.plot(time,energy,label=legende[i], c=col[i])
for i,energy2 in enumerate([kin2, pot2, kin2+pot2]):
    plt.plot(time,energy2, '--', c=col[i])
plt.legend(loc='right')
plt.xlabel('$t$ in yr')
plt.ylabel('$E$ in $10^{35}$ J')
plt.ylabel('$E$ in $10^{35}$\,J ($10^{34}$\,J dashed)')
plt.xlim(-1,80)
# Jupiter orbiting period
plt.axvline(x=11.862, color='gray', lw=0.7)
plt.axvline(x=11.862*2, color='gray', lw=0.7)
plt.axvline(x=11.862*3, color='gray', lw=0.7)
plt.axvline(x=11.862*4, color='gray', lw=0.7)
plt.axvline(x=11.862*5, color='gray', lw=0.7)
plt.axvline(x=11.862*6, color='gray', lw=0.7)
# Saturn orbiting period
plt.axvline(x=29.457, color='black', lw=0.7)
plt.axvline(x=58.914, color='black', lw=0.7)
# plt.title('Energy of the system as a function of time')
# plt.savefig(PDF/"Energy_comparison_Jupiter.pdf")


#%%
Daten_min = np.loadtxt(CSV/"sat-trajectories-to-pluto.csv",delimiter=';')
Daten_max = np.loadtxt(CSV/"sat-trajectories-to-pluto2.csv",delimiter=';')

n = int((len(Daten_min[0,:])-1)/6)    # total number of planets
number = n                        # number of planets to display

# create the variables and assign them their values via a loop
for i,name in enumerate(["x_min", "y_min", "z_min","vx_min", "vy_min", "vz_min"]):
  globals()[name] = Daten_min[:,i*n+1:(i+1)*n+1]
for i,name in enumerate(["x_max", "y_max", "z_max","vx_max", "vy_max", "vz_max"]):
  globals()[name] = Daten_max[:,i*n+1:(i+1)*n+1]

# Plot the trajectories
rcParams['figure.figsize'] = (6,3)
plt.figure()
for i in range(0,10):
    plt.plot(x_min[:,i], y_min[:,i], label=Namen[i], c=col[10-1-i])
    
plt.plot(x_min[:1850,10], y_min[:1850,10],'--', c=col[18-10])
plt.plot(x_min[:2050,11], y_min[:2050,11],'--', c=col[18-11])
plt.plot(x_min[:2650,15], y_min[:2650,15],'--', c=col[18-15])
plt.plot(x_min[:4000,16], y_min[:4000,16],'--', c=col[18-16])

for i in ((12,13,14)):
    plt.plot(x_max[:4200,i], y_max[:4200,i],'--', c=col[18-i])    
plt.xlim(-33,63)
plt.ylim(-32,47)
plt.xlabel('$x$ in AU')
plt.ylabel('$y$ in AU')
plt.legend()
plt.savefig(PDF/"different_planets.pdf")
#%%  Animate the planetary trajectories
%matplotlib auto
# col = ['#000000', '#0C5DA5', '#0C5DA5', '#0C5DA5', '#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']
col = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e', '#e377c2', '#8c564b', '#17becf', '#bcbd22']
from matplotlib import animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=120, metadata=dict(artist='Martin Beyer'), bitrate=-1)

fig = plt.figure(figsize=(6,3))
ax1 = plt.axes()
line, = ax1.plot([], [], lw=2)
plt.xlabel('$x$ in a.u.')
plt.ylabel('$y$ in a.u.')

for i in range(number):
    ax1.plot(x[:,i],y[:,i],'-', c=col[i])  

lines = []
x2 = x[:,:number]
y2 = y[:,:number]
for i in range(number):
    lobj = ax1.plot([],[],'o', label=Namen[i], c=col[i])[0]
    lines.append(lobj)
plt.legend(loc='right')
# plt.xlim(-2.1,-1.5)
# plt.ylim(-5,-4)

# plt.xlim(-2,2)
# plt.ylim(-5,2)

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
anim = animation.FuncAnimation(fig, animate, init_func=init,frames=len(x), interval=0.5, blit=True)

plt.show()
# anim.save('Test.mp4', writer=writer, dpi=400)