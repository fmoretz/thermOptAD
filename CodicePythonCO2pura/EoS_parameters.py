#-----------------------------------------------------
# This module contains constants and parameters used
# in the Equations of State for carbon dioxide
# developed by Span & Wagner and by Jäger & Span
#-----------------------------------------------------

import numpy as np
# Parameters describing the triple point of carbon dioxide (Span & Wagner, 1996)
T_t = 216.592                # In units of K
P_t = 517950                 # In units of Pa
# Parameters describing the critical point of carbon dioxide (Span & Wagner, 1996)
T_c = 304.1282               # In units of K
P_c = 7377300                # In units of Pa
r_c = 467.6                  # In units of Kg*m^(-3)
# Avogadro constant
N_a = 6.022140857E23         # In units of mol^(-1)
# Universal gas constant
R = 8.31446262             # In units of J*mol^(-1)*K^(-1)
R_m = 188.924269             # In units of J*kg^(-1)*K^(-1)
# Carbon dioxide molar mass
Mco2 = 0.0440095             # In units of Kg*mol^(-1)
# Reference state used by Jäger & Span
T_0 = 150                    # In units of K
P_0 = 101325                 # In units of Pa
# Jäger-Span parameters
g_0 = -2.6385478
g_1 = 4.5088732
g_2 = -2.0109135
g_3 = -2.7976237
g_4 = 0.26427834
g_5 = 3.8259935
g_6 = 0.31711996
g_7 = 0.0022087195
g_8 = -1.1289668
g_9 = 0.0092923982
g_10 = 3391.4617
g_0a = 0.039993365
g_1a = 0.0023945101
g_2a = 0.32839467
g_3a = 0.057918471
g_4a = 0.0023945101
g_5a = -0.0026531689
g_6a = 0.16419734
g_7a = 0.17594802
g_8a = 0.0026531689
g_0k = 0.22690751
g_1k = -0.07501975
g_2k = 0.26442913
# Span-Wagner parameters (ideal-gas part of the Helmholtz energy)
#-------------------------------------------------------
# WARNING: Span and Wagner have chosen the "ideal gas
# reference state" (h = 0 and s = 0 for ideal gas at
# T = 298.15 K and P = 101325 Pa). However, here it
# has been chosen the reference state for fluid CO2
# used in the software REFPROP (saturated liquid at
# 273.15 K, h = 200 kJ/kg, s = 1 kJ/kgK) in order to
# get the equality of the Gibbs free energy for the
# solid, liquid, and vapor phases at the triple point
#-------------------------------------------------------
a_1 = -6.12487106335325
a_2 = 5.11559631859619
a_3 = 2.5
a_4 = 1.99427042
a_5 = 0.62105248
a_6 = 0.41195293
a_7 = 1.04028922
a_8 = 0.08327678
th_4 = 3.15163
th_5 = 6.11190
th_6 = 6.77708
th_7 = 11.32384
th_8 = 27.08792
# Span-Wagner parameters (residual part of the Helmholtz energy)
n1 = np.array([0.38856823203161,2.938547594274,-5.5867188534934,-0.76753199592477,
    0.31729005580416,0.54803315897767,0.12279411220335])
d1 = np.array([1, 1, 1, 1, 2, 2, 3])
t1 = np.array([0.00, 0.75, 1.00, 2.00, 0.75, 2.00, 0.75])
n2 = np.array([2.1658961543220, 1.5841735109724,-0.23132705405503,5.8116916431436e-2,
    -0.55369137205382, 0.48946615909422,-2.4275739843501e-2,0.062494790501678,
    -0.12175860225246, -0.37055685270086,-1.6775879700426e-2,-0.11960736637987,
    -0.045619362508778, 3.5612789270346e-2,-7.4427727132052e-3,-1.7395704902432e-3,
    -0.021810121289527,0.024332166559236,-3.7440133423463e-2,0.14338715756878,
    -0.13491969083286,-0.02315122505348, 1.2363125492901e-2,0.002105832197294,
    -3.3958519026368e-4, 5.5993651771592e-3, -3.0335118055646e-4])
d2 = np.array([1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8])
t2 = np.array([1.5,1.5,2.5,0.0,1.5,2.0,0.0,1.0,2.0,3.0,6.0,3.0,6.0,8.0,6.0,0.0,7.0,12.0,
    16.0,22.0,24.0,16.0,24.0,8.0,2.0,28.0,14.0])
c2 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6])
n3 = np.array([-213.6548868832,26641.569149272,-24027.212204557,-283.41603423999,212.47284400179])
d3 = np.array([2, 2, 2, 3, 3])
t3 = np.array([1.0, 0.0, 1.0, 3.0, 3.0])
alpha3 = np.array([25, 25, 25, 15, 20])
beta3 = np.array([325, 300, 300, 275, 275])
gamma3 = np.array([1.16, 1.19, 1.19, 1.25, 1.22])
n4 = np.array([-0.66642276540751, 0.72608632349897, 0.055068668612842])
C4 = np.array([10.0, 10.0, 12.5])
alpha4 = np.array([3.5, 3.5, 3.0])
beta4 = np.array([0.3, 0.3, 1.0])
bi4 = np.array([0.875, 0.875, 0.875])
