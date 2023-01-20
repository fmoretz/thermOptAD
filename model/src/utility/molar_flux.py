import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from func import Antoine
from prop import Ant

from sklearn.preprocessing import MinMaxScaler

plt.style.use(['science', 'no-latex'])
dpi = 200
scaler = MinMaxScaler()

# evaluate the molar flux of a species ina binary mixture (A-W) where:
T = [20,  25,  30,  35, 40,  45,  50, 55]  # °C
p = [0.3, 0.4, 0.5, 0.7,  1, 1.2, 1.5, 1.6]  # atm
R = 0.0821e-3 # atm m3 / mol K
L = 5 # m first hypothesis
d = 15 # m second hypothesis

A = ['CH4', 'CO2', 'H2S']
B = 'H2O'

# safe log function
def safe_log(x):
    from math import log
    if x == 0: return   0
    else: return log(x)

# molar flux (mol/m2/d) function
def molar_flux(p, D, R, T, L, psat, phi):
    return p*D/(R*T*L) * safe_log((p-psat*phi)/(p-psat)) * (3600*24)

# water dynamic viscosity (Pas) function
def water_viscosity(T_K: float):
    A = 2.41e-5 # Pas
    B = 247.8 # K
    C = 140 # K
    return A * 10**(B/(T_K - C))

# diffusion coefficients for i species in water (m2/s)
D = {
    A[0]: [(0.0295*T[i] + 1.1246)*10**(-9) for i in range(len(T))], # Moradi et al., 2020
    A[1]: [5.35*10**(-10)*(T[i]+273.15)/((water_viscosity(T[i]+273.15)**1.035)*10)*10**(-4) for i in range(len(T))], # Rinker & Sandall (1994)
    A[2]: [1.91*10**(-9)*(T[i]+273.15)/((water_viscosity(T[i]+273.15)**0.74)*10)*10**(-4)   for i in range(len(T))]  # Halmour & Sandall (1984)
} # m2/s

# saturation pressure of species (atm)
psat = {
    A[0]: [Antoine(Ant[A[0]], T[i]+273.15) for i in range(len(T))],
    A[1]: [Antoine(Ant[A[1]], T[i]+273.15) for i in range(len(T))],
    A[2]: [Antoine(Ant[A[2]], T[i]+273.15) for i in range(len(T))],
    B:    [Antoine(Ant[B],    T[i]+273.15) for i in range(len(T))]
} # atm

# load resulting compositions (mol/mol)
df_03_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_03atm.csv', header=0)
df_04_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_04atm.csv', header=0)
df_05_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_05atm.csv', header=0)
df_07_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_07atm.csv', header=0)
df_10_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_10atm.csv', header=0)
df_12_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_12atm.csv', header=0)
df_15_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_15atm.csv', header=0)
df_16_atm = pd.read_csv(filepath_or_buffer='output/compositions/csv/xy_16atm.csv', header=0)

y_H2O = {
    '03atm': df_03_atm['yH2O'],
    '04atm': df_04_atm['yH2O'],
    '05atm': df_05_atm['yH2O'],
    '07atm': df_07_atm['yH2O'],
    '10atm': df_10_atm['yH2O'],
    '12atm': df_12_atm['yH2O'],
    '15atm': df_15_atm['yH2O'],
    '16atm': df_16_atm['yH2O'],
    }

# phi = P_H2O(T)/Psat_H2O(T) = p*y_H2O(T,p)/Psat_H2O(T)
phi = {
    '03atm': [p[0]*y_H2O['03atm'][i]/psat[B][i] for i in range(len(T))],
    '04atm': [p[1]*y_H2O['04atm'][i]/psat[B][i] for i in range(len(T))],
    '05atm': [p[2]*y_H2O['05atm'][i]/psat[B][i] for i in range(len(T))],
    '07atm': [p[3]*y_H2O['07atm'][i]/psat[B][i] for i in range(len(T))],
    '10atm': [p[4]*y_H2O['10atm'][i]/psat[B][i] for i in range(len(T))],
    '12atm': [p[5]*y_H2O['12atm'][i]/psat[B][i] for i in range(len(T))],
    '15atm': [p[6]*y_H2O['15atm'][i]/psat[B][i] for i in range(len(T))],
    '16atm': [p[7]*y_H2O['16atm'][i]/psat[B][i] for i in range(len(T))],
}

N_03atm ={
    A[0]: [molar_flux(p[0], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['03atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[0], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['03atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[0], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['03atm'][i]) for i in range(len(T))]
}
N_04atm ={
    A[0]: [molar_flux(p[1], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['04atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[1], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['04atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[1], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['04atm'][i]) for i in range(len(T))]
}
N_05atm ={
    A[0]: [molar_flux(p[2], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['05atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[2], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['05atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[2], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['05atm'][i]) for i in range(len(T))]
}
N_07atm ={
    A[0]: [molar_flux(p[3], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['07atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[3], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['07atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[3], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['07atm'][i]) for i in range(len(T))]
}
N_10atm ={
    A[0]: [molar_flux(p[4], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['10atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[4], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['10atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[4], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['10atm'][i]) for i in range(len(T))]
}
N_12atm ={
    A[0]: [molar_flux(p[5], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['12atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[5], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['12atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[5], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['12atm'][i]) for i in range(len(T))]
}
N_15atm ={
    A[0]: [molar_flux(p[6], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['15atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[6], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['15atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[6], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['15atm'][i]) for i in range(len(T))]
}
N_16atm ={
    A[0]: [molar_flux(p[7], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['16atm'][i]) for i in range(len(T))],
    A[1]: [molar_flux(p[7], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['16atm'][i]) for i in range(len(T))],
    A[2]: [molar_flux(p[7], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['16atm'][i]) for i in range(len(T))]
}

# molar fluxes evaluation (mol/m2)
S = np.pi*d**2/4

n_03atm ={
    A[0]: [molar_flux(p[0], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['03atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[0], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['03atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[0], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['03atm'][i])*S for i in range(len(T))]
}
n_04atm ={
    A[0]: [molar_flux(p[1], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['04atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[1], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['04atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[1], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['04atm'][i])*S for i in range(len(T))]
}
n_05atm ={
    A[0]: [molar_flux(p[2], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['05atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[2], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['05atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[2], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['05atm'][i])*S for i in range(len(T))]
}
n_07atm ={
    A[0]: [molar_flux(p[3], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['07atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[3], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['07atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[3], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['07atm'][i])*S for i in range(len(T))]
}
n_10atm ={
    A[0]: [molar_flux(p[4], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['10atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[4], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['10atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[4], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['10atm'][i])*S for i in range(len(T))]
}
n_12atm ={
    A[0]: [molar_flux(p[5], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['12atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[5], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['12atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[5], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['12atm'][i])*S for i in range(len(T))]
}
n_15atm ={
    A[0]: [molar_flux(p[6], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['15atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[6], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['15atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[6], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['15atm'][i])*S for i in range(len(T))]
}
n_16atm ={
    A[0]: [molar_flux(p[7], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['16atm'][i])*S for i in range(len(T))],
    A[1]: [molar_flux(p[7], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['16atm'][i])*S for i in range(len(T))],
    A[2]: [molar_flux(p[7], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['16atm'][i])*S for i in range(len(T))]
}


# visualization
ms = 2
lw = 0.75

plt.figure(dpi=dpi)
plt.subplot(2,3,1)
filename = 'absolute molar flux [mol/d]'
plt.plot(T, n_03atm[A[0]], marker='o', markersize=ms, linewidth=lw, label='P: $0.3\,atm$')
plt.plot(T, n_05atm[A[0]], marker='^', markersize=ms, linewidth=lw, label='P: $0.5\,atm$')
plt.plot(T, n_10atm[A[0]], marker='s', markersize=ms, linewidth=lw, label='P: $1.0\,atm$')
plt.plot(T, n_15atm[A[0]], marker='D', markersize=ms, linewidth=lw, label='P: $1.5\,atm$')
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
plt.ylabel(filename)
plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
plt.title('CH$_4$')

plt.subplot(2,3,2)
plt.plot(T, n_03atm[A[1]], marker='o', markersize=ms, linewidth=lw, label='P: $0.3\,atm$')
plt.plot(T, n_05atm[A[1]], marker='^', markersize=ms, linewidth=lw, label='P: $0.5\,atm$')
plt.plot(T, n_10atm[A[1]], marker='s', markersize=ms, linewidth=lw, label='P: $1.0\,atm$')
plt.plot(T, n_15atm[A[1]], marker='D', markersize=ms, linewidth=lw, label='P: $1.5\,atm$')
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
plt.title('CO$_2$')

plt.subplot(2,3,3)
plt.plot(T, n_03atm[A[2]],marker='o', markersize=ms, linewidth=lw, label='P: $0.3\,atm$')
plt.plot(T, n_05atm[A[2]],marker='^', markersize=ms, linewidth=lw, label='P: $0.5\,atm$')
plt.plot(T, n_10atm[A[2]],marker='s', markersize=ms, linewidth=lw, label='P: $1.0\,atm$')
plt.plot(T, n_15atm[A[2]],marker='D', markersize=ms, linewidth=lw, label='P: $1.5\,atm$')
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
plt.title('H$_2$S')
plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

plt.subplot(2,3,4)
filename = 'normalized molar flux [mol/d]'
plt.plot(T, n_03atm[A[0]]/np.max(n_03atm[A[0]]),marker='o', linewidth=lw, markersize=ms, label='P: $0.3\,atm$')
plt.plot(T, n_05atm[A[0]]/np.max(n_05atm[A[0]]),marker='^', linewidth=lw, markersize=ms, label='P: $0.5\,atm$')
plt.plot(T, n_10atm[A[0]]/np.max(n_10atm[A[0]]),marker='s', linewidth=lw, markersize=ms, label='P: $1.0\,atm$')
plt.plot(T, n_15atm[A[0]]/np.max(n_15atm[A[0]]),marker='D', linewidth=lw, markersize=ms, label='P: $1.5\,atm$')
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
plt.xlabel('T [°C]'); plt.ylabel(filename)
plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

plt.subplot(2,3,5)
plt.plot(T, n_03atm[A[1]]/np.max(n_03atm[A[1]]),marker='o', linewidth=lw, markersize=ms, label='P: $0.3\,atm$')
plt.plot(T, n_05atm[A[1]]/np.max(n_05atm[A[1]]),marker='^', linewidth=lw, markersize=ms, label='P: $0.5\,atm$')
plt.plot(T, n_10atm[A[1]]/np.max(n_10atm[A[1]]),marker='s', linewidth=lw, markersize=ms, label='P: $1.0\,atm$')
plt.plot(T, n_15atm[A[1]]/np.max(n_15atm[A[1]]),marker='D', linewidth=lw, markersize=ms, label='P: $1.5\,atm$')
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
plt.xlabel('T [°C]')
plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

plt.subplot(2,3,6)
plt.plot(T, n_03atm[A[2]]/np.max(n_03atm[A[2]]),marker='o', linewidth=lw, markersize=ms, label='P: $0.3\,atm$')
plt.plot(T, n_05atm[A[2]]/np.max(n_05atm[A[2]]),marker='^', linewidth=lw, markersize=ms, label='P: $0.5\,atm$')
plt.plot(T, n_10atm[A[2]]/np.max(n_10atm[A[2]]),marker='s', linewidth=lw, markersize=ms, label='P: $1.0\,atm$')
plt.plot(T, n_15atm[A[2]]/np.max(n_15atm[A[2]]),marker='D', linewidth=lw, markersize=ms, label='P: $1.5\,atm$')
plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
plt.xlabel('T [°C]')
plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
plt.savefig('molarflux.svg', dpi=dpi)
plt.show()