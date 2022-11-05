#----------------------------------------------------------
# This module contains correlation equations developed
# by Span & Wagner in order to calculate the equilibrium
# pressure for carbon dioxide at a given temperature
#----------------------------------------------------------

import numpy as np
import math as mt
# Module containing parameters used in Span-Wagner and Jäger-Span EoS
from EoS_parameters import T_t, T_c, P_t, P_c

def P_eq(m,T):
    """
    #-----------------------------------------------------
    # Function describing sublimation, melting and vapor
    # pressure (Span-Wagner's correlation equations)
    #-----------------------------------------------------
    """
    # Span-Wagner parameters (correlation equation)
    a_m = np.array([1955.5390, 2055.4593])
    a_s = np.array([-14.740846, 2.4327015, -5.3061778])
    a_v = np.array([-7.0602087, 1.9391218, -1.6463597, -3.2995634])
    t_v = np.array([1.0, 1.5, 2.0, 4.0])

    # Functions used in order to calculate the melting pressure
    diffm = T/T_t - 1
    sumM  = 1 + a_m[0]*diffm + a_m[1]*diffm**2

    # Functions used in order to calculate the sublimation pressure
    diffs = 1 - T/T_t
    sumS  = a_s[0]*diffs + a_s[1]*diffs**(1.9) + a_s[2]*diffs**(2.9)

    # Functions used in order to calculate the saturated vapor pressure
    diffv = 1 - T/T_c

    sumV = 0.
    for i in range(3):
        sumV = sumV + a_v[i]*(diffv)**t_v[i]

        #----------------------------------------
        # Set up equilibrium pressure function:
        # m = 1 --> sublimation pressure
        # m = 2 --> vapor pressure
        # m = 3 --> melting pressure
        #----------------------------------------
        if m == 1:
            p = P_t*mt.exp((T_t/T)*sumS)
        elif m == 2:
            p = P_c*mt.exp((T_c/T)*sumV)
        elif m == 3:
            p = P_t*sumM

    return p
