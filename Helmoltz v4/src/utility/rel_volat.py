# evaluation of the relative volatilities of the species
from pprint import pprint as pp
import numpy as np
import pandas as pd
from prop import kH, Ant
from func import Henry, Antoine
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import matplotlib

# matplotlib.rcParams.update({'font.size': 20, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
# matplotlib.rcParams['font.sans-serif'] = "Arial"
plt.style.use(['science', 'no-latex'])

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

dpi = 300
fig, ax = plt.subplots(dpi=dpi)
ax.plot(Temperature, rel_vol['CH4'], linewidth=1, label='CH$_4$')
ax.plot(Temperature, rel_vol['CO2'], linewidth=1, label='CO$_2$')
ax.plot(Temperature, rel_vol['H2S'], linewidth=1, label='H$_2$S')
ax.xaxis.set_ticks(np.linspace(Temperature[0], Temperature[-1], len(Temperature)+3), endpoint=True)
#et_ticks(np.linspace(min(y_H2S[len(Pressure)-1]),max(y_H2S[0]), len(Temperature)+1), endpoint=True)
ax.set_xlabel('T [°C]')
ax.set_ylabel('relative volatility [--]')
#plt.legend()
ax.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
filename = 'volatility.png'
fig.savefig(filename, dpi = dpi)

fig, ax = plt.subplots(dpi=dpi)
ax.plot(Temperature, dPsat['CH4']/np.max(dPsat['CH4']), linewidth=1, label='CH$_4$')
ax.plot(Temperature, dPsat['CO2']/np.max(dPsat['CO2']), linewidth=1, label='CO$_2$')
ax.plot(Temperature, dPsat['H2S']/np.max(dPsat['H2S']), linewidth=1, label='H$_2$S')
ax.plot(Temperature, dPsat['H2O']/np.max(dPsat['H2O']), linewidth=1, label='H$_2$O')
ax.set_xlabel('T [°C]')
ax.xaxis.set_ticks(np.linspace(Temperature[0], Temperature[-1], len(Temperature)+3), endpoint=True)
ax.set_ylabel('Rate of change [Pa/°C]')
#plt.legend()
ax.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

filename = 'rate_of_change.png'
fig.savefig(filename, dpi = dpi)

plt.show()