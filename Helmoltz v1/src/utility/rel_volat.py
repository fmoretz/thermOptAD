# evaluation of the relative volatilities of the species
from pprint import pprint as pp
import numpy as np
import pandas as pd
from prop import kH, Ant
from func import Henry, Antoine
from src.func2 import volat_y
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib

matplotlib.rcParams.update({'font.size': 20, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"

# Evaluate the relative volatility of these components at different temperature ranges
Temperature = [25, 35, 45, 55]  # °C
species = ['CH4', 'CO2', 'H2S', 'H2O']

Psat = {
    species[0]: [Antoine(Ant[species[0]], T+273.15) for T in Temperature],
    species[1]: [Antoine(Ant[species[1]], T+273.15) for T in Temperature],
    species[2]: [Antoine(Ant[species[2]], T+273.15) for T in Temperature],
    species[3]: [Antoine(Ant[species[3]], T+273.15) for T in Temperature]
}

rel_vol = {
    species[0]: [Psat[species[0]][i]/Psat[species[3]][i] for i in range(len(Temperature))],
    species[1]: [Psat[species[1]][i]/Psat[species[3]][i] for i in range(len(Temperature))],
    species[2]: [Psat[species[2]][i]/Psat[species[3]][i] for i in range(len(Temperature))],
    species[3]: [Psat[species[3]][i]/Psat[species[3]][i] for i in range(len(Temperature))]
}

# Evaluate the rate of change of vapour pressure over temperature
def dPsatdT(Psat, T):
    rate = np.empty(len(Psat))
    for i in range(1, len(Psat)):
        rate[i] = (Psat[i-1] - Psat[i])/(T[i-1] - T[i])
    return rate

dPsat = {
    species[0]: dPsatdT(Psat=Psat[species[0]], T=Temperature),
    species[1]: dPsatdT(Psat=Psat[species[1]], T=Temperature),
    species[2]: dPsatdT(Psat=Psat[species[2]], T=Temperature),
    species[3]: dPsatdT(Psat=Psat[species[3]], T=Temperature)
}

bin_vol = {
    species[0]: [volat_y(y=y[species[0]]) for i in range(len(T))],
    species[1]: [volat_y(y=y[species[1]]) for i in range(len(T))],
    species[2]: [volat_y(y=y[species[2]]) for i in range(len(T))],
    species[3]: [volat_y(y=y[species[3]]) for i in range(len(T))]
}

plt.figure()
plt.plot(Temperature, rel_vol['CH4'], linewidth=2, color='black', label='CH$_4$')
plt.plot(Temperature, rel_vol['CO2'], linewidth=2, color='red', label='CO$_2$')
plt.plot(Temperature, rel_vol['H2S'], linewidth=2, color='orange', label='H$_2$S')
plt.xlabel('Temperature [°C]')
plt.ylabel('Adimensional relative volatility ($\\alpha$) [--]')
plt.legend()
plt.grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)

plt.figure()
plt.plot(Temperature, dPsat['CH4']/np.max(dPsat['CH4']), linewidth=2, color='black', label='CH$_4$')
plt.plot(Temperature, dPsat['CO2']/np.max(dPsat['CO2']), linewidth=2, color='red', label='CO$_2$')
plt.plot(Temperature, dPsat['H2S']/np.max(dPsat['H2S']), linewidth=2, color='orange', label='H$_2$S')
plt.plot(Temperature, dPsat['H2O']/np.max(dPsat['H2O']), linewidth=2, color='blue', label='H$_2$O')
plt.xlabel('Temperature [°C]')
plt.ylabel('Rate of change ($dP_{sat}/dT$) [Pa/°C]')
plt.legend()
plt.grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)

plt.show()