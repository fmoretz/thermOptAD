import sys
import numpy as np
import pylab
from pylab import *
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import matplotlib.ticker as ticker
import pandas as pd
import openpyxl
from openpyxl import load_workbook

# Reading data set (if you have more than one sheet, change the name of the sheet yuo want analyse )
path = 'Data.xlsx'
data  = pd.read_excel(path, sheet_name="P_sub")
data1 = pd.read_excel(path, sheet_name="P_melting")
data2 = pd.read_excel(path, sheet_name="P_sat")
# sorting data set
df =data.sort_values("T [K]")
df1 = data1.sort_values("T [K]")
df2 = data2.sort_values("T [K]")

# Set up the features of plot 
plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
matplotlib.rcParams.update({'font.size': 12, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['font.sans-serif'] = "Arial"

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Set x and y values
x  = df.get('T [K]')
x1 = df1.get('T [K]')
x2 = df2.get('T [K]')
y1 = df.get('P [MPa]')
y2 = df.get('P_eq [MPa]')
y3 = df1.get('P [MPa]')
y4 = df1.get('P_eq [MPa]')
y5 = df2.get('P [MPa]')
y6 = df2.get('P_eq [MPa]')

axes[0].scatter(x, y1, s=35, c='b', alpha=0.5, lw=2.0, edgecolors='b', label='Exp. data')
axes[0].plot(x, y2,'b', lw=2.0,label='Sublimation line (model)')
axes[0].scatter(x1, y3, s=35, c='g', alpha=0.5, lw=2.0, edgecolors='g', label='Exp. data')
axes[0].plot(x1, y4,'g',lw=2.0,label='Melting line (model)')
axes[0].scatter(x2, y5, s=35, c='r', alpha=0.5, lw=2.0, edgecolors='r', label='Exp. data')
axes[0].plot(x2, y6,'r',lw=2.0,label='Saturation line (model)')

axes[0].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
axes[0].set_xlabel(r"$T\,[\mathrm{K}]$", fontsize=15)
axes[0].set_ylabel(r"$P\,[\mathrm{MPa}]$", fontsize=15)
axes[0].set_xlim(170, 305)
axes[0].set_ylim(0, 270)
axes[0].legend(loc='upper left', fontsize=12)

#-------------------------------------------------------------------------------------------------------------------
# Percentage deviation diagram
#-------------------------------------------------------------------------------------------------------------------
axes[1].scatter(x, 100*(y2-y1)/y1, s=35, c='b', alpha=0.5, lw=2.0, edgecolors='b', label='Sub. pressure deviation')
axes[1].scatter(x1, 100*(y4-y3)/y3, s=35, c='g', alpha=0.5, lw=2.0, edgecolors='g', label='Melt. pressure deviation')
axes[1].scatter(x2, 100*(y6-y5)/y5, s=35, c='r', alpha=0.5, lw=2.0, edgecolors='r', label='Sat. pressure deviation')
axes[1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
axes[1].axhline(y=0, c='k', linestyle='-', lw=1.2)
axes[1].set_xlabel(r"$T\,[\mathrm{K}]$", fontsize=15)
axes[1].set_ylabel(r"$E_{r}\,[\mathrm{\%}]$", fontsize=15)
axes[1].set_xlim(165, 305)
axes[1].set_ylim(-0.6, 0.6)
axes[1].legend(loc='lower right', fontsize=12)

fig.savefig('Result.png',dpi=200)
plt.show()