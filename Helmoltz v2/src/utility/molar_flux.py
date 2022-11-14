import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from func import Antoine
from prop import Ant
from pprint import pprint as pp

plt.style.use(['science', 'ieee', 'no-latex'])
dpi = 300

# evaluate the molar flux of a species ina binary mixture (A-W) where:
T = [20,  25,  30,  35, 40,  45,  50, 55]  # °C
p = [0.3, 0.4, 0.5, 0.7,  1, 1.2, 1.5, 1.6]  # atm
R = 0.0821e-3 # atm m3 / mol K
L = 5 # m first hypothesis

A = ['CH4', 'CO2', 'H2S']
B = 'H2O'

# safe log function
def safe_log(x):
    from math import log
    if x == 0: return   0
    else: return log(x)

# molar flux (mol/m3/s) function
def molar_flux(p, D, R, T, L, psat, phi):
    return p*D/(R*T*L) * safe_log((p-psat*phi)/(p-psat))  

# water dynamic viscosity (Pas) function
def water_viscosity(T_K: float):
    """Water viscosity calculation

    Args:
        T_K (float): temperature of the water in K

    Returns:
        float: visocosity of the water in Pas
    """
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

# phi = {
#     A[0]:,
#     A[1]:,
#     A[2]:    
# }

# N ={
#     A[0]: [molar_flux(p[i], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi[i]) for i in range(len(T))],
#     A[1]: [molar_flux(p[i], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi[i]) for i in range(len(T))],
#     A[2]: [molar_flux(p[i], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi[i]) for i in range(len(T))]
# } # mol/m2/s