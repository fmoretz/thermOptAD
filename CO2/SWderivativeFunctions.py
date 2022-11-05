#------------------------------------------------
# This module contains the Helmholtz energy and
# its derivatives, which are necessary for the
# calculation of phase equilibrium conditions
# and of properties from Span & Wagner's model
#------------------------------------------------

import math as mt
import numpy as np
from scipy.misc import derivative
# Module containing parameters used in Span-Wagner and Jäger-Span EoS
from EoS_parameters import *

# Span-Wagner functions
a0  = np.array([a_4, a_5, a_6, a_7, a_8])
th0 = np.array([th_4, th_5, th_6, th_7, th_8])

def Phi(tau,d_r):
    """
    #--------------------------------------------------------
    # Function describing the dimensionless Helmholtz energy
    #--------------------------------------------------------
    """

    # Ideal-gas part of the Helmholtz energy
    B = mt.log(d_r) + a_1 + a_2*tau + a_3*mt.log(tau)
    B0 = 0.
    for i in range(4):
        B0 = B0 + a0[i]*mt.log(1-mt.exp(-tau*th0[i]))

    # Calculation of the first term of the Helmholtz energy's residual part
    sum1 = 0.
    for i in range(6):
        sum1 = sum1 + n1[i]*(d_r**d1[i])*(tau**t1[i])

    # Calculation of the second term of the Helmholtz energy's residual part
    sum2 = 0.
    for i in range(26):
        sum2 = sum2 + n2[i]*(d_r**d2[i])*(tau**t2[i])*exp(-d_r**c2[i])

    # Calculation of the third term of the Helmholtz energy's residual part
    sum3 = 0.
    for i in range(4):
        sum3 = sum3 + n3[i]*(d_r**d3[i])*(tau**t3[i])*exp(-alpha3[i]*(d_r-1)**2-beta3[i]*(tau-gamma3[i])**2)

    # Calculation of the fourth term of the Helmholtz energy's residual part
    sum4 = 0.
    for i in range(2):
        delta = (1-tau + 0.7*((d_r-1)**2)**(1/0.6))**2 + beta4[i]*((d_r-1)**2)**alpha4[i]
        sum4 = sum4 + n4[i]*(delta**bi4[i])*d_r*exp(-(C4[i]*(d_r-1)**2)-275*((tau-1)**2))

    # Dimensionless Helmholtz energy
    return (B + B0 + sum1 + sum2 + sum3 + sum4)

def dphi0_t(tau):
    """
    #----------------------------------------------
    # Function describing the derivative of the
    # ideal-gas part of the Helmholtz energy with
    # respect to inverse reduced temperature
    #----------------------------------------------
    """

    sum1 = 0.
    for i in range(4):
        sum1 = sum1 + a0[i]*th0[i]*(1/(1-mt.exp(-tau*th0[i])) - 1)

    # Derivative of the Ideal-gas part of the Helmholtz energy with respect to tau
    return (a_2 + a_3/tau + sum1)

def dphi0_tt(tau):
    """
    #-----------------------------------------------
    # Function describing the derivative of dphi0_t
    # with respect to inverse reduced temperature
    #-----------------------------------------------
    """

    # Derivative of dphi0_t with respect to tau
    return derivative(dphi0_t, tau)

def dphir_t(tau,d_r):
    """
    #---------------------------------------------
    # Function describing the derivative of the
    # residual part of the Helmholtz energy with
    # respect to inverse reduced temperature
    #---------------------------------------------
    """

    diff = d_r - 1

    # Derivative of the first term of Helmholtz energy's residual part with respect to tau
    sum1 = 0.
    for i in range(6):
        sum1 = sum1 + n1[i]*t1[i]*(d_r**d1[i])*(tau**(t1[i] - 1))

    # Derivative of the second term of Helmholtz energy's residual part with respect to tau
    sum2 = 0.
    for j in range(26):
        sum2 = sum2 + n2[j]*t2[j]*(d_r**d2[j])*(tau**(t2[j] -1))*mt.exp(-d_r**c2[j])

    # Derivative of the third term of Helmholtz energy's residual part with respect to tau
    sum3 = 0.
    for z in range(4):
        sum3 = sum3 + n3[z]*(d_r**d3[z])*(tau**t3[z])*mt.exp(-(alpha3[z]*diff**2)\
        - beta3[z]*((tau-gamma3[z])**2))*((t3[z]/tau) - 2*beta3[z]*(tau-gamma3[z]))

    # Derivative of the fourth term of Helmholtz energy's residual part with respect to tau
    sum4 = 0.
    for h in range(2):
        theta    = 1 - tau + 0.7*(diff**2)**(1/0.6)
        delta    = theta**2 + beta4[h]*(diff**2)**alpha4[h]
        deltab_t = -2*bi4[h]*theta*(delta**(bi4[h]-1))
        psi      = mt.exp(-(C4[h]*diff**2)-275*((tau - 1)**2))
        psi_t    = -550*(tau - 1)*psi
        sum4     = sum4 + n4[h]*d_r*((delta**bi4[h])*psi_t + psi*deltab_t)

    # Derivative of the residual part of the Helmholtz energy with respect to tau
    return (sum1 + sum2 + sum3 + sum4)

def dphir_tt(tau,d_r):
    """
    #-----------------------------------------------
    # Function describing the derivative of dphir_t
    # with respect to inverse reduced temperature
    #-----------------------------------------------
    """

    diff = d_r - 1

    # Derivative of the first term of dphir_t with respect to tau
    sum1 = 0.
    for i in range(6):
        sum1 = sum1 + n1[i]*t1[i]*(t1[i] - 1)*(d_r**d1[i])*(tau**(t1[i] - 2))

    # Derivative of the second term of dphir_t with respect to tau
    sum2 = 0.
    for j in range(26):
        sum2 = sum2 + n2[j]*t2[j]*(t2[j]-1)*(d_r**d2[j])*(tau**(t2[j]-2))*mt.exp(-d_r**c2[j])

    # Derivative of the third term of dphir_t with respect to tau
    sum3 = 0.
    for z in range(4):
        eTau = tau-gamma3[z]
        sum3 = sum3 + n3[z]*(d_r**d3[z])*(tau**t3[z])*mt.exp(-(alpha3[z]*diff**2)\
        -beta3[z]*eTau**2)*((t3[z]/tau-2*beta3[z]*eTau)**2-t3[z]/(tau**2)-2*beta3[z])

    # Derivative of the fourth term of dphir_t with respect to tau
    sum4 = 0.
    for h in range(2):
        theta     = 1 - tau + 0.7*(diff**2)**(1/0.6)
        delta     = theta**2 + beta4[h]*(diff**2)**alpha4[h]
        deltab_t  = -2*bi4[h]*theta*(delta**(bi4[h]-1))
        deltab_tt = 2*bi4[h]*(delta**(bi4[h]-1))+4*(delta**(bi4[h]-2))*bi4[h]*(bi4[h]-1)*theta**2
        psi       = mt.exp(-(C4[h]*diff**2) - 275*((tau - 1)**2))
        psi_t     = -550*(tau - 1)*psi
        psi_tt    = 550*(-550*(tau - 1)**2 - 1)*psi
        sum4      = sum4 + n4[h]*d_r*((delta**bi4[h])*psi_tt + 2*psi_t*deltab_t + psi*deltab_tt)

    # Derivative of dphir_t energy with respect to tau
    return (sum1 + sum2 + sum3 + sum4)

def dphi_r(d_r,tau):
    """
    #-----------------------------------------
    # Function describing the derivative of 
    # the residual part of the Helmholtz
	# energy with respect to reduced density
    #------------------------------------------
    """

    diff = d_r - 1

    # Derivative of the first term of Helmholtz energy's residual part with respect to d_r
    sum1 = 0.
    for i in range(6):
        sum1 = sum1 + n1[i]*d1[i]*(d_r**(d1[i]-1))*(tau**t1[i])

    # Derivative of the second term of Helmholtz energy's residual part with respect to d_r
    sum2 = 0.
    for j in range(26):
        sum2 = sum2 + n2[j]*mt.exp(-d_r**c2[j])*(d_r**(d2[j]-1))*(tau**t2[j])*(d2[j]-c2[j]*d_r**c2[j])

    # Derivative of the third term of Helmholtz energy's residual part with respect to d_r
    sum3 = 0.
    for z in range(4):
        sum3 = sum3 + n3[z]*(d_r**d3[z])*(tau**t3[z])*mt.exp(-(alpha3[z]*diff**2)\
        - beta3[z]*((tau-gamma3[z])**2))*((d3[z]/d_r) - 2*alpha3[z]*diff)

    # Derivative of the fourth term of Helmholtz energy's residual part with respect to d_r
    sum4 = 0.
    for h in range(2):
        theta  = 1 - tau + 0.7*(diff**2)**(1/0.6)
        delta  = theta**2 + beta4[h]*(diff**2)**alpha4[h]
        delta1 = diff*((1.4*theta/0.3)*(diff**2)**(1/0.6-1) + 2*beta4[h]*alpha4[h]*(diff**2)**(alpha4[h]-1))
        deltab = bi4[h]*(delta**(bi4[h]-1)*delta1)
        psi    = mt.exp(-(C4[h]*diff**2) - 275*((tau-1)**2))
        psi1   = -2*C4[h]*diff*psi
        sum4   = sum4 + n4[h]*((delta**bi4[h])*(psi + d_r*psi1) + d_r*psi*deltab)

    # Derivative of the residual part of the Helmholtz energy with respect to reduced density
    return (sum1 + sum2 + sum3 + sum4)

def d2phi_r(d_r,tau):
    """
    #-----------------------------------------
    # Function describing the derivative of 
    # dphi_r with respect to reduced density
    #-----------------------------------------
    """

    diff = d_r - 1

    # Derivative of the first term of dphi_r with respect to reduced density
    sum1 = 0.
    for i in range(6):
        sum1 = sum1 + n1[i]*d1[i]*(d1[i]-1)*(d_r**(d1[i]-2))*(tau**t1[i])

    # Derivative of the second term of dphi_r with respect to reduced density
    sum2 = 0.
    for j in range(26): 
       dr   = d_r**c2[j]
       sum2 = sum2 + n2[j]*mt.exp(-dr)*(d_r**(d2[j]-2))*(tau**t2[j])*((d2[j]-c2[j]*dr)*(d2[j]-1-c2[j]*dr)-dr*c2[j]**2)

    # Derivative of the third term of dphi_r with respect to reduced density
    sum3 = 0.
    for z in range(4):
        eTau = tau-gamma3[z]
        dr2  = d_r**d3[z]
        d3z  = d3[z] - 1
        sum3 = sum3 + n3[z]*(tau**t3[z])*mt.exp(-(alpha3[z]*diff**2)-beta3[z]*eTau**2)*(-2*alpha3[z]*dr2\
            + 4*(alpha3[z]**2)*dr2*diff**2 -4*d3[z]*alpha3[z]*diff*(dr2/d_r) + d3[z]*d3z*d_r**(d3z-1))

    # Derivative of the fourth term of dphi_r with respect to reduced density
    sum4 = 0.
    for h in range(2):
        theta   = 1 - tau + 0.7*(diff**2)**(1/0.6)
        delta   = theta**2 + beta4[h]*(diff**2)**alpha4[h]
        d1_1    = (diff**2)**(1/0.6 - 1)
        delta1  = diff*((1.4*theta/0.3)*d1_1 + 2*beta4[h]*alpha4[h]*(diff**2)**(alpha4[h]-1))
        d2_1    = 4*beta4[h]*alpha4[h]*(alpha4[h]-1)*(diff**2)**(alpha4[h]-2)
        d2_2    = (2.8*theta/0.3)*(1/0.6 - 1)*(diff**2)**(1/0.6-2)
        delta2  = (delta1/diff) + (diff**2)*(d2_1 + (0.98/0.09)*d1_1**2 + d2_2)
        deltab  = bi4[h]*(delta**(bi4[h]-1)*delta1)
        deltab2 = bi4[h]*(delta**(bi4[h]-1)*delta2 + (bi4[h]-1)*(delta**(bi4[h]-2))*delta1**2)
        psi     = mt.exp(-(C4[h]*diff**2) - 275*((tau-1)**2))
        psi1    = -2*C4[h]*diff*psi
        psi2    = 2*C4[h]*psi*(2*C4[h]*diff**2 - 1)
        sum4    = sum4 + n4[h]*((delta**bi4[h])*(2*psi1 + d_r*psi2)+2*deltab*(psi + d_r*psi1)+d_r*psi*deltab2)

    return (sum1 + sum2 + sum3 + sum4)

def dphir_rt(tau,d_r):
    """
    #----------------------------------------------
    # Function describing the derivative of dphi_r 
    # with respect to inverse reduced temperature
    #---------------------------------------------
    """

    diff = d_r - 1

    # Derivative of the first term of dphi_r with respect to tau
    sum1 = 0.
    for i in range(6):
        sum1 = sum1 + n1[i]*d1[i]*t1[i]*(d_r**(d1[i]-1))*(tau**(t1[i] - 1))

    # Derivative of the second term of dphi_r with respect to tau
    sum2 = 0.
    for j in range(26):
        sum2 = sum2 + n2[j]*t2[j]*mt.exp(-d_r**c2[j])*(d_r**(d2[j]-1))*(tau**(t2[j]-1))*(d2[j]-c2[j]*d_r**c2[j])

    # Derivative of the third term of dphi_r with respect to tau
    sum3 = 0.
    for z in range(4):
        eTau = tau-gamma3[z]
        sum3  = sum3 + n3[z]*(d_r**d3[z])*(tau**t3[z])*mt.exp(-(alpha3[z]*diff**2)\
           -beta3[z]*eTau**2)*((d3[z]/d_r)-2*alpha3[z]*diff)*((t3[z]/tau)-2*beta3[z]*eTau)

    # Derivative of the fourth term of dphi_r with respect to tau
    sum4 = 0.
    for h in range(2):
        dbi4     = delta**bi4[h]
        theta    = 1 - tau + 0.7*(diff**2)**(1/0.6)
        delta    = theta**2 + beta4[h]*(diff**2)**alpha4[h]
        delta1   = diff*((1.4*theta/0.3)*(diff**2)**(1/0.6-1)+2*beta4[h]*alpha4[h]*(diff**2)**(alpha4[h]-1))
        deltab   = bi4[h]*(delta**(bi4[h]-1)*delta1)
        deltab_t = -2*bi4[h]*theta*(delta**(bi4[h]-1))
        delta_t  = -2*bi4[h]*theta*(bi4[h]-1)*(delta**(bi4[h]-2))*delta1
        delta_rt = -(delta**(bi4[h]-1))*(1.4*bi4[h]/0.3)*diff*(diff**2)**(1/0.6-1) + delta_t
        psi      = mt.exp(-(C4[h]*diff**2) - 275*((tau-1)**2))
        psi1     = -2*C4[h]*diff*psi
        psi_t    = -550*(tau -1)*psi
        psi_rt   = 1100*C4(h)*diff*(tau-1)*psi
        sum4     = sum4 + n4[h]*(dbi4*(psi_t + d_r*psi_rt) + d_r*psi_t*deltab + deltab_t*(psi+d_r*psi1)+d_r*psi*delta_rt)

    # Derivative of dphi_r with respect to tau
    return (sum1 + sum2 + sum3 + sum4)
