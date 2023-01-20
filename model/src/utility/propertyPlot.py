import sys
import numpy as np
import pylab
from pylab import *
from func import *
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker

def pressure_densityPlot(specie):

  #importing necessary properties
  from prop import R_m, T_c, r_c
  from func import dphi_r
  T_c = T_c[specie]
  r_c = r_c[specie]
  R_m = R_m[specie]

  # Set up the features of plot
  plt.style.use(['science', 'no-latex'])

  dpi = 500

  fig, axes = plt.subplots(1, 2, dpi=dpi, figsize=(10, 5))

  T_range = [25+273, 30+273, 35+273, 40+273, 45+273, 50+273, 55+273] # K

  tau  = [T_c/T_range[0], T_c/T_range[1], T_c/T_range[2], T_c/T_range[3], T_c/T_range[4], T_c/T_range[5], T_c/T_range[6]]

  # Set x and y values
  x  = np.arange(1200)
  y1 = np.empty(shape=(len(x),))
  y2 = np.empty(shape=(len(x),))
  y3 = np.empty(shape=(len(x),))
  y4 = np.empty(shape=(len(x),))
  y5 = np.empty(shape=(len(x),))
  y6 = np.empty(shape=(len(x),))
  y7 = np.empty(shape=(len(x),))

  for i in range(1200):
    y1[i] = x[i]*R_m*T_range[0]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[0], specie))
    y2[i] = x[i]*R_m*T_range[1]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[1], specie))
    y3[i] = x[i]*R_m*T_range[2]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[2], specie))
    y4[i] = x[i]*R_m*T_range[3]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[3], specie))
    y5[i] = x[i]*R_m*T_range[4]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[4], specie))
    y6[i] = x[i]*R_m*T_range[5]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[5], specie))
    y7[i] = x[i]*R_m*T_range[6]*(1 + (x[i]/r_c)*dphi_r(x[i]/r_c, tau[5], specie))

  labels = [
    'T = ' + str(T_range[0]) + ' K',
    'T = ' + str(T_range[1]) + ' K',
    'T = ' + str(T_range[2]) + ' K',
    'T = ' + str(T_range[3]) + ' K',
    'T = ' + str(T_range[4]) + ' K',
    'T = ' + str(T_range[5]) + ' K',
    'T = ' + str(T_range[6]) + ' K',
  ]

  axes[0].plot(x, y1/10**6,'k',     label=labels[0])
  axes[0].plot(x, y2/10**6,'b',     label=labels[1])
  axes[0].plot(x, y3/10**6,'c',     label=labels[2])
  axes[0].plot(x, y4/10**6,'m',     label=labels[3])
  axes[0].plot(x, y5/10**6,'r',     label=labels[4])
  axes[0].plot(x, y6/10**6,'orange',label=labels[5])
  axes[0].plot(x, y7/10**6,'lime',  label=labels[6])

  axes[0].grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
  axes[0].set_xlabel(r"$\mathrm{\rho_{H_2S}[kg/m^3]}$")
  axes[0].set_ylabel(r"$\mathrm{P\,[MPa]}$")
  axes[0].set_xlim(0, 1000)
  axes[0].xaxis.set_ticks(np.linspace(0, 1000, 11, endpoint=True))
  axes[0].set_ylim(0, 50)
  axes[0].yaxis.set_ticks(np.linspace(0, 50, 11, endpoint=True))
  #axes[0].legend(loc='upper left')

  axes[1].plot(x, y1/10**6,'k',     label=labels[0])
  axes[1].plot(x, y2/10**6,'b',     label=labels[1])
  axes[1].plot(x, y3/10**6,'c',     label=labels[2])
  axes[1].plot(x, y4/10**6,'m',     label=labels[3])
  axes[1].plot(x, y5/10**6,'r',     label=labels[4])
  axes[1].plot(x, y6/10**6,'orange',label=labels[5])
  axes[1].plot(x, y7/10**6,'lime',  label=labels[6])

  axes[1].grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
  axes[1].set_xlabel(r"$\mathrm{\rho_{H_2O}\,[kg/m^3]}$")
  axes[1].set_ylabel(r"$\mathrm{P\,[MPa]}$")
  axes[1].set_xlim(800, 860)
  axes[1].xaxis.set_ticks(np.linspace(800, 860, 11, endpoint=True))
  axes[1].set_ylim(0, 10)
  axes[1].yaxis.set_ticks(np.linspace(0, 10, 11,  endpoint=True))
  #axes[1].legend(loc='upper left')

  filename1 = str(specie) + '_densityPlot.svg'
  filename2 = str(specie) + '_densityPlot.png'
  fig.savefig(filename1, dpi = dpi)
  fig.savefig(filename2, dpi = dpi)
  return print('\nfigure saved succesfully!\n')

# ==================================================================== #
# def temperature_fugacityPlot(specie):

#   #importing necessary properties
#   from prop import R_m, T_c, r_c, d_vap, d_liq
#   from func import dphi_r
#   T_c = T_c[specie]
#   r_c = r_c[specie]
#   R_m = R_m[specie]

#   # Set up the features of plot
#   plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
#   matplotlib.rcParams.update({'font.size': 12, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
#   matplotlib.rcParams['font.sans-serif'] = "Arial"

#   fig, axes = plt.subplots(1, 1, figsize=(20, 10))

#   P  = [0.1, 0.2, 0.5, 0.7, 1.0, 1.2, 1.5] # atm

#   # Set x and y values
#   x  = np.arange(300, 351) # Temperature gradient
#   tau = np.empty(shape=(len(x),))
#   y1 = np.empty(shape=(len(x),))
#   y2 = np.empty(shape=(len(x),))
#   y3 = np.empty(shape=(len(x),))
#   y4 = np.empty(shape=(len(x),))
#   y5 = np.empty(shape=(len(x),))
#   y6 = np.empty(shape=(len(x),))
#   y7 = np.empty(shape=(len(x),))

#   for i in range(len(x)):

#     d_r1 = density(d_vap[specie], x[i], P[0]*101325, specie)/r_c
#     d_r2 = density(d_vap[specie], x[i], P[1]*101325, specie)/r_c
#     d_r3 = density(d_vap[specie], x[i], P[2]*101325, specie)/r_c
#     d_r4 = density(d_vap[specie], x[i], P[3]*101325, specie)/r_c
#     d_r5 = density(d_vap[specie], x[i], P[4]*101325, specie)/r_c
#     d_r6 = density(d_vap[specie], x[i], P[5]*101325, specie)/r_c
#     d_r7 = density(d_vap[specie], x[i], P[6]*101325, specie)/r_c

#     tau[i] = T_c/x[i]

#     y1[i] = np.exp(phiR(tau[i],d_r1,specie) +d_r1*dphi_r(d_r1,tau[i],specie)-np.log(1 + d_r1*dphi_r(d_r1,tau[i],specie)))
#     y2[i] = np.exp(phiR(tau[i],d_r2,specie) +d_r2*dphi_r(d_r2,tau[i],specie)-np.log(1 + d_r2*dphi_r(d_r2,tau[i],specie)))
#     y3[i] = np.exp(phiR(tau[i],d_r3,specie) +d_r3*dphi_r(d_r3,tau[i],specie)-np.log(1 + d_r3*dphi_r(d_r3,tau[i],specie)))
#     y4[i] = np.exp(phiR(tau[i],d_r4,specie) +d_r4*dphi_r(d_r4,tau[i],specie)-np.log(1 + d_r4*dphi_r(d_r4,tau[i],specie)))
#     y5[i] = np.exp(phiR(tau[i],d_r5,specie) +d_r5*dphi_r(d_r5,tau[i],specie)-np.log(1 + d_r5*dphi_r(d_r5,tau[i],specie)))
#     y6[i] = np.exp(phiR(tau[i],d_r6,specie) +d_r6*dphi_r(d_r6,tau[i],specie)-np.log(1 + d_r6*dphi_r(d_r6,tau[i],specie)))
#     y7[i] = np.exp(phiR(tau[i],d_r7,specie) +d_r7*dphi_r(d_r7,tau[i],specie)-np.log(1 + d_r7*dphi_r(d_r7,tau[i],specie)))

#   labels = [
#     'P = ' + str(P[0]) + ' atm',
#     'P = ' + str(P[1]) + ' atm',
#     'P = ' + str(P[2]) + ' atm',
#     'P = ' + str(P[3]) + ' atm',
#     'P = ' + str(P[4]) + ' atm',
#     'P = ' + str(P[5]) + ' atm',
#     'P = ' + str(P[6]) + ' atm'
#   ]

#   axes.plot(x, y1,'k',     lw=2.0,label=labels[0])
#   axes.plot(x, y2,'b',     lw=2.0,label=labels[1])
#   axes.plot(x, y3,'c',     lw=2.0,label=labels[2])
#   axes.plot(x, y4,'m',     lw=2.0,label=labels[3])
#   axes.plot(x, y5,'r',     lw=2.0,label=labels[4])
#   axes.plot(x, y6,'orange',lw=2.0,label=labels[5])
#   axes.plot(x, y7,'lime',  lw=2.0,label=labels[6])

#   axes.grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
#   axes.set_xlabel(r"$T\,[\mathrm{K}]$", fontsize=15)
#   axes.set_ylabel(r"$\phi_{CH_4}[\mathrm{--}]$", fontsize=15)
#   lim_min = min(x); lim_max = 250
#   axes.set_xlim(lim_min, lim_max)
#   axes.xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
#   lim_min = float(min(y1))/1.004; lim_max = float(max(y7))*1.003
#   axes.set_ylim(lim_min, lim_max)
#   axes.yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
#   #axes.legend(loc='upper left', fontsize=12)

#   # axes[1].plot(x, y1,'k',     lw=2.0,label=labels[0])
#   # axes[1].plot(x, y2,'b',     lw=2.0,label=labels[1])
#   # axes[1].plot(x, y3,'c',     lw=2.0,label=labels[2])
#   # axes[1].plot(x, y4,'m',     lw=2.0,label=labels[3])
#   # axes[1].plot(x, y5,'r',     lw=2.0,label=labels[4])
#   # axes[1].plot(x, y6,'orange',lw=2.0,label=labels[5])
#   # axes[1].plot(x, y7,'lime',  lw=2.0,label=labels[6])

#   # axes[1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
#   # axes[1].set_xlabel(r"$T\,[\mathrm{K}]$", fontsize=15)
#   # axes[1].set_ylabel(r"$\phi\,[\mathrm{--}]$", fontsize=15)
#   # # axes[1].set_xlim(0, 500)
#   # # axes[1].xaxis.set_ticks(np.linspace(0, 500, 11, endpoint=True))
#   # # axes[1].set_ylim(0, 10)
#   # # axes[1].yaxis.set_ticks(np.linspace(0, 10, 11, endpoint=True))
#   # axes[1].legend(loc='upper left', fontsize=12)

#   filename = str(specie) + '_fugacityPlot.svg'
#   fig.savefig(filename, dpi=800)
#   return print('\nfigure saved succesfully!\n')

#! Remember to change xlim-ylim ranges of the plots
species     = ['CH4', 'CO2', 'H2S', 'H2O']
pressure_densityPlot(species[2])