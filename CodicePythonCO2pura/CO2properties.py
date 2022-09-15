#-------------------------------------------
# This module contains functions which are
# necessary for the calculation of carbon
# dioxide thermo-chemical properties at a
# given pressure and temperature
#-------------------------------------------

import numpy as np
import math as mt
# Module containing the derivatives of the Helmholtz energy
import SWderivativeFunctions as SWdF

def VLProperties(q,T,d):
    """
    #--------------------------------------------
    # Function describing  liquid/vapor phase
    # properties of carbon dioxide
    #--------------------------------------------
    """

    # Set up reduced density and the inverse reduced temperature
    d_r = d/r_c
    tau = T_c/T

    # Functions used in order to calculate SW's heat capacity
    dpdv = 1 + 2*d_r*SWdF.dphi_r(d_r,tau) + (d_r**2)*SWdF.d2phi_r(d_r,tau)
    dpdt = 1 + d_r*SWdF.dphi_r(d_r,tau) - tau*d_r*SWdF.dphir_rt(tau,d_r)

    # Dimensionless SW's Helmholtz energy
    Ar      = SWdF.Phi(tau,d_r)
    # Molar entropy
    entropy = R*(tau*(SWdF.dphi0_t(tau) + SWdF.dphir_t(tau,d_r)) - SWdF.Phi(tau,d_r))
    # Molar enthalpy
    enthalpy = R*T*(1+tau*(SWdF.dphi0_t(tau)+SWdF.dphir_t(tau,d_r))+d_r*SWdF.dphi_r(d_r,tau))
    # Isobaric heat capacity
    cp = -R*((tau**2)*(SWdF.dphi0_tt(tau) + SWdF.dphir_tt(tau,d_r)) - (dpdt**2)/dpdv)
    # Isochoric heat capacity
    cv = -R*(tau**2)*(SWdF.dphi0_tt(tau) + SWdF.dphir_tt(tau,d_r))
    # Cubic expansion coefficient
    cubicExpCoefficient = dpdt/(dpdv*T)
    # Surface tension
    surfTension = 0.07863*((1-T/T_c)**1.254)

    if q == 0
        return Ar
    elif q == 1
        return entropy
    elif q == 2
        return enthalpy
    elif q == 3
        return cp
    elif q == 4
        return cv
    elif q == 5
        return cubicExpCoefficient
    elif q == 6
        return ThermalConductivity(tau,d_r)
    elif q == 7
        return Viscosity(T,d)
    elif q == 8
        return surfTension

def SolidProperties(q,T,P):
    """
    #---------------------------------------------
    # Function describing  solid phase properties
    # of carbon dioxide (Jäger-Span's equations)
    #---------------------------------------------
    """

    # Functions used in order to calculate JS's molar Gibbs free energy and its derivatives
    theta = T/T_0
    Pi    = P/P_0
    f_a1 = g_1a*mt.log((theta**2 - g_2a*theta + g_3a)/(1 - g_2a + g_3a))
    f_a2 = g_4a*mt.log((theta**2 + g_2a*theta + g_3a)/(1 + g_2a + g_3a))
    f_a3 = g_5a*(mt.atan((theta - g_6a)/g_7a) - mt.atan((1 - g_6a)/g_7a))
    f_a4 = g_8a*(mt.atan((theta + g_6a)/g_7a) - mt.atan((1 + g_6a)/g_7a))
    K_t  = g_0k*theta**2 + g_1k*theta + g_2k
    f_a  = g_0a*(theta**2 - 1) + f_a1 + f_a2 +f_a3 + f_a4

    df_a1 = g_1a*((2*theta - g_2a)/(theta**2 - g_2a*theta + g_3a))
    df_a2 = g_4a*((2*theta + g_2a)/(theta**2 + g_2a*theta + g_3a))
    df_a3 = (g_5a/g_7a)*(1/(1 + ((theta - g_6a)/g_7a)**2))
    df_a4 = (g_8a/g_7a)*(1/(1 + ((theta + g_6a)/g_7a)**2))
    # Derivative of f_a with respect to temperature
    df_a = 2*g_0a*theta + df_a1 + df_a2 + df_a3 + df_a4

    d2f_a1 = g_1a*((2*(theta**2 - g_2a*theta + g_3a)-(2*theta - g_2a)**2)/(theta**2 - g_2a*theta + g_3a)**2)
    d2f_a2 = g_4a*((2*(theta**2 + g_2a*theta + g_3a)-(2*theta + g_2a)**2)/(theta**2 + g_2a*theta + g_3a)**2)
    d2f_a3 = - (g_5a/g_7a**2)*(2/(1 + ((theta - g_6a)/g_7a)**2)**2)*((theta - g_6a)/g_7a)
    d2f_a4 = - (g_8a/g_7a**2)*(2/(1 + ((theta + g_6a)/g_7a)**2)**2)*((theta + g_6a)/g_7a)
    # Derivative of df_a with respect to temperature
    d2f_a = 2*g_0a + d2f_a1 + d2f_a2 + d2f_a3 + d2f_a4

    # Functions used in order to calculate JS's molar Gibbs free energy
    G1 = g_3*(mt.log((theta**2 + g_4**2)/(1+g_4**2))-(2*theta/g_4)*(mt.atan(theta/g_4)-mt.atan(1/g_4)))
    G2 = g_5*(mt.log((theta**2 + g_6**2)/(1+g_6**2))-(2*theta/g_6)*(mt.atan(theta/g_6)-mt.atan(1/g_6)))
    G3 = g_7*(Pi-1)*(mt.exp(f_a) + K_t*g_8) + g_9*K_t*((Pi + g_10)**(6/7) - (1 + g_10)**(6/7))

    # Functions used in order to calculate JS's molar entropy
    E = 2*(g_3/g_4)*(mt.atan(theta/g_4)-mt.atan(1/g_4)) + 2*(g_5/g_6)*(mt.atan(theta/g_6)-mt.atan(1/g_6))
    E1 = g_7*(Pi-1)*(mt.exp(f_a)*df_a + g_8*(2*g_0k*theta + g_1k))
    E2 = g_9*(2*g_0k*theta + g_1k)*((Pi + g_10)**(6/7) - (1 + g_10)**(6/7))
    E1cec = g_7*(mt.exp(f_a)*df_a + g_8*(2*g_0k*theta + g_1k))

    # Functions used in order to calculate JS's heat capacity
    Cp1 = 2*g_3/(g_4**2 + theta**2) + 2*g_5/(g_6**2 + theta**2)
    Cp2 = g_7*(Pi-1)*(mt.exp(f_a)*d2f_a + mt.exp(f_a)*df_a**2 + g_8*2*g_0k)
    Cp3 = g_9*2*g_0k*((Pi + g_10)**(6/7) - (1 + g_10)**(6/7))

    # JS's molar Gibbs free energy
    g = R*T_0*(g_0 + g_1*(theta - 1) + g_2*(theta - 1)**2 + G1 +G2 + G3)
    # Molar volume computed using the derivative with respect to pressure of JS's Gibbs free energy
    mVolume = ((R*T_0)/P_0)*(g_7*(mt.exp(f_a) + g_8*K_t) + g_9*K_t*(6/7)*((Pi + g_10)**(-1/7)))
    # Molar entropy computed using the derivative with respect to temperature of JS's Gibbs free energy
    entropy = - R*(g_1 + 2*g_2*(theta - 1) - E + E1 + E2)
    # Molar enthalpy computed using Gibbs free energy and entropy
    enthalpy = g + T*entropy
    # Isobaric heat capacity computed using the derivative with respect to temperature of entropy
    cp = - ((R*T)/T_0)*(2*g_2 - Cp1 + Cp2 + Cp3)
    # Isochoric heat capacity computed using the derivative with respect to temperature of internal energy
    cv = cp - P*(R/P_0)*(E1cec + g_9*(2*g_0k*theta + g_1k)*(6/7)*((Pi + g_10)**(-1/7)))
    # Cubic expansion coefficient computed using molar volume and its derivative with respect to temperature
    cubicExpCoefficient = (R/P_0)*(E1cec + g_9*(2*g_0k*theta + g_1k)*(6/7)*((Pi + g_10)**(-1/7)))/mVolume

    if q == 0
        return g
    elif q == 1
        return Mco2/mVolume
    elif q == 2
        return entropy
    elif q == 3
        return enthalpy
    elif q == 4
        return cp
    elif q == 5
        return cv
    elif q == 6
        return cubicExpCoefficient

def ThermalConductivity(tau,d_r):
    """
    #-------------------------------------------
    # Function describing  liquid/vapor phase
    # thermal conductivity of carbon dioxide
    #-------------------------------------------
    """

    l0 = np.array([0.028067404, 0.022856419, -0.0074162421])
    B1 = np.array([0.0100128, 0.0560488, -0.081162, 0.0624337, -0.0206336, 0.00253248])
    B2 = np.array([0.00430829, -0.0358563, 0.067148, -0.0522855, 0.0174571, -0.00196414])

    # Set up reduced temperature
    t_r = 1/tau

    # Functions used in order to calculate Huber's thermal conductivity
    dT = t_r -1
    dD = d_r -1

    suml0 = 0.0151874307
    for i in range(2):
        suml0 = suml0 + l0[i]/(t_r**i)

    # First term of Huber's thermal conductivity
    lambda0 = mt.sqrt(t_r)/suml0

    # Second term of Huber's thermal conductivity
    lambda1 = 0.
    for j in range(5): 
        lambda1 = lambda1 + (B1[j]+B2[j]*t_r)*(d_r**j)

    # Third term of Huber's thermal conductivity
    lambda2 = (-17.47-44.88*dT)/(0.8563-mt.exp(8.865*dT+4.16*dD**2 +2.302*dT*dD -dD**3)-0.4503*dD-7.197*dT)

    # Thermal conductivity
    return (lambda1 + (lambda0 + lambda2)/1000)

def Viscosity(T,d):
    """
    #------------------------------------
    # Function describing  liquid/vapor 
    # phase viscosity of carbon dioxide
    #------------------------------------
    """

    V0 = 1749.35489318835
    V1 = -369.06930000712
    V2 = 5423856.34887691
    V3 = -2.2128385216835
    V4 =-269503.247933569
    V5 = 73145.0215318260
    V6 = 5.34368649509278
    n = np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.5, 5.5])
    B = np.array([219.73999,-1015.3226,2471.0125,-3375.1717,2491.6597,-787.26086,14.085455,-0.34664158])

    # Set up reduced temperatures and reduced density
    t_s = T/200.76
    t_r = T/T_t
    d_r = d/1178.53

    # Functions used in order to calculate Laesecke's viscosity
    etaFactor = 1000*(1178.53**(2/3))*mt.sqrt(R*T_t)/((Mco2**(1/6))*N_a**(1/3))

    sumB = -19.572881
    for i in range(7):
        sumB = sumB + B[i]/(t_s**n[i])

    # First term of Laesecke's viscosity
    eta0 = 1.0055*mt.sqrt(T)/(V0+V1*T**(1/6)+V2*mt.exp(V3*T**(1/3))+(V4+V5*T**(1/3))/mt.exp(T**(1/3))+V6*mt.sqrt(T))

    # Second term of Laesecke's viscosity
    eta1 = (eta0*sumB*(0.378421E-09)**3)*N_a/Mco2

    # Third term of Laesecke's viscosity
    eta2 = etaFactor*(0.360603235428487*t_r*d_r**3 + (d_r**2 + d_r**8.06282737481277)/(t_r-0.121550806591497))

    # Thermal conductivity
    return ((eta0 + d*eta1 + eta2)/1000)

