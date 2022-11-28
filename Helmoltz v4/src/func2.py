import mpmath as mp
import math as mt
import warnings
import numpy as np
from scipy import optimize as opt
from scipy.optimize import fsolve
from prop2 import *
import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path
plt.style.use(['science', 'no-latex'])
warnings.filterwarnings('ignore')

dps = 5


def phi0(tau, d_r, species):
    """
    #--------------------------------------------------
    # Function describing the Ideal-gas part of
    # the dimensionless Helmholtz energy
    # d_r = reduced density (rho/rho_crit)
    # tau = inverse reduced temperature (Tcrit/T)
    # rho in Kg*m^(-3)
    # T in K
    #---------------------------------------------------
    """

    from prop2 import a0, th0, a, T_c, r_c, OX, R_m

    if species == 'O2':
        d0 = 101325 / (r_c[species] * R_m[species] * 298.15)
        phi_id = safe_log(d_r / d0) + OX[0] * tau ** 1.5 + OX[1] / tau ** 2 + OX[2] * safe_log(tau) \
                 + OX[3] * tau + OX[4] * safe_log(mp.exp(tau * OX[6]) - 1) + OX[5] \
                 * safe_log(2 / 3 * mp.exp(-tau * OX[7]) + 1) + OX[8]
    else:
        B = a[species][4] * tau ** (1 / 3) + a[species][5] / tau ** (3 / 2) + a[species][6] / tau ** (7 / 4)
        B0 = 0.
        for i in range(len(a0[species])):
            B0 = B0 + a0[species][i] * mp.log(1 - mp.exp(-tau * th0[species][i]))

        phi_id = safe_log(d_r) + a[species][0] + a[species][1] * tau + a[species][2] * safe_log(tau) + a[species][3] * (
                T_c[species] / tau) ** 1.5 + B0 + B

    # Ideal-gas part of the dimensionless Helmholtz energy
    return phi_id


def phiR(tau, d_r, species):
    """
    #--------------------------------------------------------
    # Function describing the dimensionless Helmholtz
    # energy's residual part
    # d_r = reduced density (rho/rho_crit)
    # tau = inverse reduced temperature (Tcrit/T)
    # rho in Kg*m^(-3)
    # T in K
    #--------------------------------------------------------
    """

    from prop2 import n1, n2, n3, n4
    from prop2 import d1, t1, d2, t2, c2, d3, t3, alpha3, beta3, gamma3, DD3
    from prop2 import Ai4, alpha4, beta4, bi4, C4, Di4

    # Calculation of the first term of the Helmholtz energy's residual part
    sum1 = 0.
    for i in range(len(n1[species])):
        sum1 = sum1 + n1[species][i] * (d_r ** d1[species][i]) * (tau ** t1[species][i])

    # Calculation of the second term of the Helmholtz energy's residual part
    sum2 = 0.
    for i in range(len(n2[species])):
        sum2 = sum2 + n2[species][i] * (d_r ** d2[species][i]) * (tau ** t2[species][i]) * mt.exp(
            -d_r ** c2[species][i])

    # Calculation of the third term of the Helmholtz energy's residual part
    sum3 = 0.
    for i in range(len(n3[species])):
        sum3 = sum3 + n3[species][i] * (d_r ** d3[species][i]) * (tau ** t3[species][i]) * mt.exp(
            -alpha3[species][i] * (d_r - DD3[species][i]) ** 2 - beta3[species][i] * (tau - gamma3[species][i]) ** 2)

    # Calculation of the fourth term of the Helmholtz energy's residual part
    sum4 = 0.
    for i in range(len(n4[species])):
        delta = (1 - tau + Ai4[species][i] * ((d_r - 1) ** 2) ** (1 / 0.6)) ** 2 + beta4[species][i] * (
                (d_r - 1) ** 2) ** alpha4[species][i]
        sum4 = sum4 + n4[species][i] * (delta ** bi4[species][i]) * d_r * mt.exp(
            -(C4[species][i] * (d_r - 1) ** 2) - Di4[species][i] * ((tau - 1) ** 2))

        # Dimensionless Helmholtz energy
    return (sum1 + sum2 + sum3 + sum4)


def Phi(tau, d_r, species):
    """
    #--------------------------------------------------------
    # Function describing the dimensionless Helmholtz energy
    # d_r = reduced density (rho/rho_crit)
    # tau = inverse reduced temperature (Tcrit/T)
    # rho in Kg*m^(-3)
    # T in K
    #--------------------------------------------------------
    """

    # Ideal and residual part of the Helmholtz energy
    Id_term = phi0(tau, d_r, species)
    Res_term = phiR(tau, d_r, species)

    # Dimensionless Helmholtz energy
    return (Id_term + Res_term)


def dphi0_t(tau, species):
    """
    #----------------------------------------------
    # Function describing the derivative of the
    # ideal-gas part of the Helmholtz energy with
    # respect to inverse reduced temperature
    #----------------------------------------------
    """

    from prop2 import th0, a0, a, T_c, OX

    if species == 'O2':
        dp_id = 1.5 * OX[0] * tau ** 0.5 - 2 * OX[1] / tau ** 3 + OX[2] / tau + OX[
            4] * OX[6] / (1 - mp.exp(-tau * OX[6])) + OX[3] - (2 / 3) * OX[5] * \
                OX[7] * mp.exp(-tau * OX[7]) / (2 / 3 * mp.exp(-tau * OX[7]) + 1)
    else:
        B = 1 / 3 * a[species][4] / tau ** (2 / 3) - 1.5 * a[species][5] / tau ** (2.5) - 1.75 * a[species][
            6] / tau ** (2.75)
        sum1 = 0.
        for i in range(len(a0[species])):
            sum1 = sum1 + a0[species][i] * th0[species][i] * (1 / (1 - mt.exp(-tau * th0[species][i])) - 1)

        dp_id = a[species][1] + a[species][2] / tau - (
                1.5 * a[species][3] * T_c[species] ** 1.5) / tau ** 2.5 + sum1 + B

    # Derivative of the Ideal-gas part of the Helmholtz energy with respect to tau
    return dp_id


def dphi0_tt(tau, species):
    """
    #-----------------------------------------------
    # Function describing the derivative of dphi0_t
    # with respect to inverse reduced temperature
    #-----------------------------------------------
    """

    from prop2 import a0, th0, a, T_c, OX

    if species == 'O2':
        d2p_id = 0.75 * OX[0] / tau ** 0.5 - OX[2] / tau ** 2 + OX[4] * (
                OX[6] ** 2) * mp.exp(tau * OX[6]) / (mp.exp(tau * OX[6]) - 1) + 6 * \
                 OX[1] / tau ** 4 + (2 / 3) * OX[5] * (OX[7] ** 2) * mp.exp(
            -tau * OX[7]) / (2 / 3 * mp.exp(-tau * OX[7]) + 1) ** 2
    else:
        B = -2 / 9 * a[species][4] / tau ** (5 / 3) + 3.75 * a[species][5] / tau ** (3.5) + 4.8125 * a[species][
            6] / tau ** (3.75)
        sum1 = 0.
        for i in range(len(a0[species])):
            sum1 = sum1 + a0[species][i] * (th0[species][i] ** 2) * mt.exp(-tau * th0[species][i]) * (
                    (1 - mt.exp(-tau * th0[species][i])) ** (-2))

        d2p_id = B - a[species][2] / (tau ** 2) - sum1 + 2.5 * (1.5 * a[species][3] * T_c[species] ** 1.5) / tau ** 3.5

    # Derivative of dphi0_t with respect to tau
    return d2p_id


def dphir_t(tau, d_r, species):
    """
    #---------------------------------------------
    # Function describing the derivative of the
    # residual part of the Helmholtz energy with
    # respect to inverse reduced temperature
    #---------------------------------------------
    """

    diff = d_r - 1

    from prop2 import n1, n2, n3, n4
    from prop2 import d1, t1, d2, t2, c2, d3, t3, alpha3, beta3, gamma3, DD3
    from prop2 import Ai4, alpha4, beta4, bi4, C4, Di4

    # Derivative of the first term of Helmholtz energy's residual part with respect to tau
    sum1 = 0.
    for i in range(len(n1[species])):
        sum1 = sum1 + n1[species][i] * t1[species][i] * (d_r ** d1[species][i]) * (tau ** (t1[species][i] - 1))

    # Derivative of the second term of Helmholtz energy's residual part with respect to tau
    sum2 = 0.
    for j in range(len(n2[species])):
        sum2 = sum2 + n2[species][j] * t2[species][j] * (d_r ** d2[species][j]) * (
                tau ** (t2[species][j] - 1)) * mt.exp(-d_r ** c2[species][j])

    # Derivative of the third term of Helmholtz energy's residual part with respect to tau
    sum3 = 0.
    for z in range(len(n3[species])):
        sum3 = sum3 + n3[species][z] * (d_r ** d3[species][z]) * (tau ** t3[species][z]) * mt.exp(
            -(alpha3[species][z] * (d_r - DD3[species][z]) ** 2) - beta3[species][z] * (
                    (tau - gamma3[species][z]) ** 2)) * (
                       (t3[species][z] / tau) - 2 * beta3[species][z] * (tau - gamma3[species][z]))

    # Derivative of the fourth term of Helmholtz energy's residual part with respect to tau
    sum4 = 0.
    for h in range(len(n4[species])):
        theta = 1 - tau + Ai4[species][h] * (diff ** 2) ** (1 / 0.6)
        delta = theta ** 2 + beta4[species][h] * (diff ** 2) ** alpha4[species][h]
        deltab_t = -2 * bi4[species][h] * theta * (delta ** (bi4[species][h] - 1))
        psi = mt.exp(-(C4[species][h] * diff ** 2) - Di4[species][h] * ((tau - 1) ** 2))
        psi_t = -2 * Di4[species][h] * (tau - 1) * psi
        sum4 = sum4 + n4[species][h] * d_r * ((delta ** bi4[species][h]) * psi_t + psi * deltab_t)

    # Derivative of the residual part of the Helmholtz energy with respect to tau
    return (sum1 + sum2 + sum3 + sum4)


def dphir_tt(tau, d_r, species):
    """
    #-----------------------------------------------
    # Function describing the derivative of dphir_t
    # with respect to inverse reduced temperature
    #-----------------------------------------------
    """

    from prop2 import n1, n2, n3
    from prop2 import d1, t1, d2, t2, c2, d3, t3, alpha3, beta3, gamma3, DD3
    from prop2 import Ai4, alpha4, beta4, bi4, C4, Di4

    # Derivative of the first term of dphir_t with respect to tau
    sum1 = 0.
    for i in range(len(n1[species])):
        sum1 = sum1 + n1[species][i] * t1[species][i] * (t1[species][i] - 1) * (d_r ** d1[species][i]) * (
                tau ** (t1[species][i] - 2))

    # Derivative of the second term of dphir_t with respect to tau
    sum2 = 0.
    for j in range(len(n2[species])):
        sum2 = sum2 + n2[species][j] * t2[species][j] * (t2[species][j] - 1) * (d_r ** d2[species][j]) * (
                tau ** (t2[species][j] - 2)) * mp.exp(-d_r ** c2[species][j])

    # Derivative of the third term of dphir_t with respect to tau
    sum3 = 0.
    for z in range(len(n3[species])):
        diff = d_r - DD3[species][z]
        eTau = tau - gamma3[species][z]
        sum3 = sum3 + n3[species][z] * (d_r ** d3[species][z]) * (tau ** t3[species][z]) * mt.exp(
             -(alpha3[species][z] * diff ** 2) - beta3[species][z] * eTau ** 2) * (
                   (t3[species][z] / tau - 2 * beta3[species][z] * eTau) ** 2 - t3[species][z] / (tau ** 2) - 2 *
                   beta3[species][z])

    # Derivative of the fourth term of dphir_t with respect to tau
    sum4 = 0.
    for h in range(len(n4[species])):
        diff  = d_r - 1
        theta = 1 - tau + Ai4[species][h] * (diff ** 2) ** (1 / 0.6)
        delta = theta ** 2 + beta4[species][h] * (diff ** 2) ** alpha4[species][h]
        deltab_t = -2 * bi4[species][h] * theta * (delta ** (bi4[species][h] - 1))
        deltab_tt = 2 * bi4[species][h] * (delta ** (bi4[species][h] - 1)) + 4 * (delta ** (bi4[species][h] - 2)) * \
                   bi4[species][h] * (bi4[species][h] - 1) * theta ** 2
        psi = mt.exp(-(C4[species][h] * diff ** 2) - Di4[species][h] * ((tau - 1) ** 2))
        psi_t = -2 * Di4[species][h] * (tau - 1) * psi
        psi_tt = 2 * Di4[species][h] * (2 * Di4[species][h] * (tau - 1) ** 2 - 1) * psi
        sum4 = sum4 + n4[species][h] * d_r * ((delta ** bi4[species][h]) * psi_tt + 2 * psi_t * deltab_t + psi * deltab_tt)

    # Derivative of dphir_t energy with respect to tau
    return (sum1 + sum2 + sum3 + sum4)


def dphi_r(d_r, tau, species):
    """
    #-----------------------------------------
    # Function describing the derivative of
    # the residual part of the Helmholtz
    # energy with respect to reduced density
    #------------------------------------------
    """

    from prop2 import n1, n2, n3, n4
    from prop2 import d1, t1, d2, t2, c2, d3, t3, alpha3, beta3, gamma3, DD3
    from prop2 import Ai4, alpha4, beta4, bi4, C4, Di4

    # Derivative of the first term of Helmholtz energy's residual part with respect to d_r
    sum1 = 0.
    for i in range(len(n1[species])):
        sum1 = sum1 + n1[species][i] * d1[species][i] * (d_r ** (d1[species][i] - 1)) * (tau ** t1[species][i])

    # Derivative of the second term of Helmholtz energy's residual part with respect to d_r
    sum2 = 0.
    for j in range(len(n2[species])):
        sum2 = sum2 + n2[species][j] * mt.exp(-d_r ** c2[species][j]) * (d_r ** (d2[species][j] - 1)) * (
                tau ** t2[species][j]) * (d2[species][j] - c2[species][j] * d_r ** c2[species][j])

    # Derivative of the third term of Helmholtz energy's residual part with respect to d_r
    sum3 = 0.
    for z in range(len(n3[species])):
        diff = d_r - DD3[species][z]
        sum3 = sum3 + n3[species][z] * (d_r ** d3[species][z]) * (tau ** t3[species][z]) * mt.exp(
        -(alpha3[species][z] * diff ** 2) - beta3[species][z] * ((tau - gamma3[species][z]) ** 2)) * (
                   (d3[species][z] / d_r) - 2 * alpha3[species][z] * diff)

    # Derivative of the fourth term of Helmholtz energy's residual part with respect to d_r
    sum4 = 0.
    for h in range(len(n4[species])):
        diff = d_r - 1
    theta = 1 - tau + Ai4[species][h] * (diff ** 2) ** (1 / 0.6)
    delta = theta ** 2 + beta4[species][h] * (diff ** 2) ** alpha4[species][h]
    delta1 = diff * ((2 * Ai4[species][h] * theta / 0.3) * (diff ** 2) ** (1 / 0.6 - 1) + 2 * beta4[species][h] *
                     alpha4[species][h] * (diff ** 2) ** (alpha4[species][h] - 1))
    deltab = bi4[species][h] * (delta ** (bi4[species][h] - 1) * delta1)
    psi = mt.exp(-(C4[species][h] * diff ** 2) - Di4[species][h] * ((tau - 1) ** 2))
    psi1 = -2 * C4[species][h] * diff * psi
    sum4 = sum4 + n4[species][h] * ((delta ** bi4[species][h]) * (psi + d_r * psi1) + d_r * psi * deltab)

    # Derivative of the residual part of the Helmholtz energy with respect to reduced density
    return (sum1 + sum2 + sum3 + sum4)


def d2phi_r(d_r, tau, species):
    """
    #-----------------------------------------
    # Function describing the derivative of
    # dphi_r with respect to reduced density
    #-----------------------------------------
    """

    from prop2 import n1, n2, n3
    from prop2 import d1, t1, d2, t2, c2, d3, t3, alpha3, beta3, gamma3, DD3
    from prop2 import Ai4, alpha4, beta4, bi4, C4, Di4

    # Derivative of the first term of dphi_r with respect to reduced density
    sum1 = 0.
    for i in range(len(n1[species])):
        sum1 = sum1 + n1[species][i] * d1[species][i] * (d1[species][i] - 1) * (d_r ** (d1[species][i] - 2)) * (
                tau ** t1[species][i])

    # Derivative of the second term of dphi_r with respect to reduced density
    sum2 = 0.
    for j in range(len(n2[species])):
        dr = d_r ** c2[species][j]
        sum2 = sum2 + n2[species][j] * mt.exp(-dr) * (d_r ** (d2[species][j] - 2)) * (tau ** t2[species][j]) * (
                (d2[species][j] - c2[species][j] * dr) * (d2[species][j] - 1 - c2[species][j] * dr) - dr * c2[species][
            j] ** 2)

    # Derivative of the third term of dphi_r with respect to reduced density
    sum3 = 0.
    for z in range(len(n3[species])):
        diff = d_r - DD3[species][z]
        eTau = tau - gamma3[species][z]
        dr2 = d_r ** d3[species][z]
        d3z = d3[species][z] - 1
        sum3 = sum3 + n3[species][z] * (tau ** t3[species][z]) * mt.exp(
            -(alpha3[species][z] * diff ** 2) - beta3[species][z] * eTau ** 2) * (
                       -2 * alpha3[species][z] * dr2 + 4 * (alpha3[species][z] ** 2) * dr2 * diff ** 2 - 4 *
                       d3[species][z] * alpha3[species][z] * diff * (dr2 / d_r) + d3[species][z] * d3z * d_r ** (
                               d3z - 1))

    # Derivative of the fourth term of dphi_r with respect to reduced density
    sum4 = 0.
    for h in range(len(n4[species])):
        diff = d_r - 1
        theta = 1 - tau + Ai4[species][h] * (diff ** 2) ** (1 / 0.6)
        delta = theta ** 2 + beta4[species][h] * (diff ** 2) ** alpha4[species][h]
        d1_1 = (diff ** 2) ** (1 / 0.6 - 1)
        delta1 = diff * (
                (2 * Ai4[species][h] * theta / 0.3) * d1_1 + 2 * beta4[species][h] * alpha4[species][h] * (
                    diff ** 2) ** (
                        alpha4[species][h] - 1))
        d2_1 = 4 * beta4[species][h] * alpha4[species][h] * (alpha4[species][h] - 1) * (diff ** 2) ** (
                    alpha4[species][h] - 2)
        d2_2 = (4 * Ai4[species][h] * theta / 0.3) * (1 / 0.6 - 1) * (diff ** 2) ** (1 / 0.6 - 2)
        delta2 = (delta1 / diff) + (diff ** 2) * (d2_1 + (2 * (Ai4[species][h] ** 2) / 0.09) * d1_1 ** 2 + d2_2)
        deltab = bi4[species][h] * (delta ** (bi4[species][h] - 1) * delta1)
        deltab2 = bi4[species][h] * (delta ** (bi4[species][h] - 1) * delta2 + (bi4[species][h] - 1) * (
                delta ** (bi4[species][h] - 2)) * delta1 ** 2)
        psi = mt.exp(-(C4[species][h] * diff ** 2) - Di4[species][h] * ((tau - 1) ** 2))
        psi1 = -2 * C4[species][h] * diff * psi
        psi2 = 2 * C4[species][h] * psi * (2 * C4[species][h] * diff ** 2 - 1)
        sum4 = sum4 + n4[species][h] * ((delta ** bi4[species][h]) * (2 * psi1 + d_r * psi2) + 2 * deltab * (
                psi + d_r * psi1) + d_r * psi * deltab2)

    return (sum1 + sum2 + sum3 + sum4)


def dphir_rt(tau, d_r, species):
    """
    #----------------------------------------------
    # Function describing the derivative of dphi_r
    # with respect to inverse reduced temperature
    #---------------------------------------------
    """

    from prop2 import n1, n2, n3
    from prop2 import d1, t1, d2, t2, c2, d3, t3, alpha3, beta3, gamma3, DD3
    from prop2 import Ai4, alpha4, beta4, bi4, C4, Di4

    # Derivative of the first term of dphi_r with respect to tau
    sum1 = 0.
    for i in range(len(n1[species])):
        sum1 = sum1 + n1[species][i] * d1[species][i] * t1[species][i] * (d_r ** (d1[species][i] - 1)) * (
                tau ** (t1[species][i] - 1))

    # Derivative of the second term of dphi_r with respect to tau
    sum2 = 0.
    for j in range(len(n2[species])):
        sum2 = sum2 + n2[species][j] * t2[species][j] * mt.exp(-d_r ** c2[species][j]) * (
                d_r ** (d2[species][j] - 1)) * (tau ** (t2[species][j] - 1)) * (
                       d2[species][j] - c2[species][j] * d_r ** c2[species][j])

    # Derivative of the third term of dphi_r with respect to tau
    sum3 = 0.
    for z in range(len(n3[species])):
        diff = d_r - DD3[species][z]
        eTau = tau - gamma3[species][z]
        sum3 = sum3 + n3[species][z] * (d_r ** d3[species][z]) * (tau ** t3[species][z]) * mt.exp(
              -(alpha3[species][z] * diff ** 2) - beta3[species][z] * eTau ** 2) * (
                   (d3[species][z] / d_r) - 2 * alpha3[species][z] * diff) * (
                   (t3[species][z] / tau) - 2 * beta3[species][z] * eTau)

    # Derivative of the fourth term of dphi_r with respect to tau
    sum4 = 0.
    for h in range(len(n4[species])):
        diff = d_r - 1
        theta = 1 - tau + Ai4[species][h] * (diff ** 2) ** (1 / 0.6)
        delta = theta ** 2 + beta4[species][h] * (diff ** 2) ** alpha4[species][h]
        dbi4 = delta ** bi4[species][h]
        delta1 = diff * ((2 * Ai4[species][h] * theta / 0.3) * (diff ** 2) ** (1 / 0.6 - 1) + 2 * beta4[species][h] *
                         alpha4[species][h] * (diff ** 2) ** (alpha4[species][h] - 1))
        deltab = bi4[species][h] * (delta ** (bi4[species][h] - 1) * delta1)
        deltab_t = -2 * bi4[species][h] * theta * (delta ** (bi4[species][h] - 1))
        delta_t = -2 * bi4[species][h] * theta * (bi4[species][h] - 1) * (delta ** (bi4[species][h] - 2)) * delta1
        delta_rt = -(delta ** (bi4[species][h] - 1)) * (2 * Ai4[species][h] * bi4[species][h] / 0.3) * diff * (
                    diff ** 2) ** (
                           1 / 0.6 - 1) + delta_t
        psi = mt.exp(-(C4[species][h] * diff ** 2) - Di4[species][h] * ((tau - 1) ** 2))
        psi1 = -2 * C4[species][h] * diff * psi
        psi_t = -2 * Di4[species][h] * (tau - 1) * psi
        psi_rt = 4 * Di4[species][h] * C4[species][h] * diff * (tau - 1) * psi
        sum4 = sum4 + n4[species][h] * (
                dbi4 * (psi_t + d_r * psi_rt) + d_r * psi_t * deltab + deltab_t * (
                    psi + d_r * psi1) + d_r * psi * delta_rt)

    # Derivative of dphi_r with respect to tau
    return (sum1 + sum2 + sum3 + sum4)


def VLProperties(q, T, d, species):
    """
    #--------------------------------------------
    # Function describing  liquid/vapor phase
    # properties of a compound
    #--------------------------------------------
    """

    from prop2 import r_c
    from prop2 import T_c
    from prop2 import R

    # Set up reduced density and the inverse reduced temperature
    d_r = d / r_c[species]
    tau = T_c[species] / T

    # Functions used in order to calculate SW's heat capacity
    dpdv = 1 + 2 * d_r * dphi_r(d_r, tau, species) + (d_r ** 2) * d2phi_r(d_r, tau, species)
    dpdt = 1 + d_r * dphi_r(d_r, tau, species) - tau * d_r * dphir_rt(tau, d_r, species)

    # Dimensionless SW's Helmholtz energy
    Ar = Phi(tau, d_r, species)
    # Molar entropy
    entropy = R * (tau * (dphi0_t(tau, species) + dphir_t(tau, d_r, species)) - Phi(tau, d_r, species))
    # Molar enthalpy
    enthalpy = R * T * (
            1 + tau * (dphi0_t(tau, species) + dphir_t(tau, d_r, species)) + d_r * dphi_r(d_r, tau, species))
    # Isobaric heat capacity
    cp = -R * ((tau ** 2) * (dphi0_tt(tau, species) + dphir_tt(tau, d_r, species)) - (dpdt ** 2) / dpdv)
    # Isochoric heat capacity
    cv = -R * (tau ** 2) * (dphi0_tt(tau, species) + dphir_tt(tau, d_r, species))
    # Cubic expansion coefficient
    cubicExpCoefficient = dpdt / (dpdv * T)
    # Fugacity coeffcient
    phi_coef = mt.exp(
        phiR(tau, d_r, species) + d_r * dphi_r(d_r, tau, species) - np.log(1 + d_r * dphi_r(d_r, tau, species)))

    if q == 0:
        return Ar
    elif q == 1:
        return entropy
    elif q == 2:
        return enthalpy
    elif q == 3:
        return cp
    elif q == 4:
        return cv
    elif q == 5:
        return cubicExpCoefficient
    elif q == 6:
        return phi_coef
# ------------------------------------------------------------------------------------------------------


def density(x0, T, P, species):
    """
    #----------------------------------------------
    # Function providing density of carbon dioxide
    # at a given pressure and temperature
    # T in K
    # P in Pa
    #----------------------------------------------
    """

    # Set up start value of density in units of Kg*m^(-3)
    # if m==1:
    #    x0 = d_liq
    # elif m==2:
    #    x0 = d_vap

    # ---------------------------------------------
    # Newton-Raphson algorithm used to calculate
    # the zero of a non-analytic function
    # ---------------------------------------------

    x = opt.newton(rho, x0, fprime=drho, args=(T, P, species), tol=1.0E-6, maxiter=int(1e6))

    return x


def rho(x, T, P, species):
    """
    #---------------------------------------------------
    # Function describing relation between pressure and
    # density: P = RT*rho*(1 + rho_ridotta*(phi_r)')
    # T in K
    # P in Pa
    #---------------------------------------------------
    """
    from prop2 import r_c
    from prop2 import T_c
    from prop2 import R_m

    # Set up reduced density and the inverse reduced temperature
    d_r = x / r_c[species]
    tau = T_c[species] / T

    # Derivative of the residual part of the Helmholtz energy with respect to reduced density
    dphi = dphi_r(d_r, tau, species)

    # -----------------------------------------------------
    # Function rho(x) --> variable x = density [Kg/m^(3)]
    # -----------------------------------------------------
    return x * R_m[species] * T * (1 + d_r * dphi) - P


def drho(x, T, P, species):
    """
    #-----------------------------------------------
    # Function describing the derivative of rho(x)
    # T in K
    # P in Pa
    #-----------------------------------------------
    """
    from prop2 import r_c
    from prop2 import T_c
    from prop2 import R_m

    # Set up reduced density and the inverse reduced temperature
    d_r = x / r_c[species]
    tau = T_c[species] / T

    # Derivatives of the residual part of the Helmholtz energy with respect to reduced density
    dphi = dphi_r(d_r, tau, species)
    ddphi_r = d2phi_r(d_r, tau, species)

    # ------------------------------------------------------
    # Function rho'(x) --> variable x = density [Kg/m^(3)]
    # ------------------------------------------------------
    return R_m[species] * T * (1 + 2 * d_r * dphi + (d_r ** 2) * ddphi_r)


def Antoine(x, T):
    '''
    Normal Antoine Equation
    log10Pev = A - B/(T+C)
    A, B and C = fitted coefficients
    Pvap = saturation pressure in bar
    T = temperature in K
    return: pressure of the species in atm
    '''
    return (10 ** (x[0] - x[1] / (T + x[2]))) * 0.986923


def Henry(x, T):
    '''
    Extended Henry's Equation
    lnH = A + B/T + ClnT + D*T
    A, B, C, and D = fitted coefficients
    H = henry's constants in kPa
    T = temperature in K
    return: henry's constants of the species in kPa
    '''
    return np.exp(x[0] + x[1] / T + x[2] * np.log(T) + x[3] * T)


def RR(alpha, z, k):
    '''
    Rachford-Rice equation, solution
    for ideal gas/ideal liquid flash systems
    z: feed molar composition
    k: pressure ratio, defined as Psat/Ptot
    Psat: saturation pressure of the species
    Ptot: total pressure of the system
    return: after zero this function, it returns
    the value of alpha (vapor fraction)
    '''
    zdk = np.empty(shape=(len(z),))
    dk = np.empty(shape=(len(z),))

    for i in range(len(z)):
        dk[i] = (k[i] - 1)
        zdk[i] = z[i] * dk[i]

    return sum(zdk / (1 + alpha * dk))


def root_Zed(w, Tcr, Pcr, P, T, Psat, display=False, u=2, v=-1):
    '''
    Find the compressibility factor
    for a species with Peng-Robinson EOS.
    objective function:
    0 = Zed^3 + alpha*Zed^2 + beta*Zed + gamma
    where:
    alpha = - 1 - B + u*B
    beta  = A + v*B^2 - u*B - u*B^2
    gamma = -A*B - v*B^2 - v*B^3
    S = 0.37464 + 1.5422 * w - 0.26992 * w**2
    k = ( 1 + S*(1-(T/Tcr)**0.5) )**2
    a = (0.45724*k*(Rgas*Tcr)**2)/Pcr
    b = (0.07780*Rgas*Tcr)/Pcr
    A = (a*P)/(Rgas*T)**2
    B = (b*P)/(Rgas*T)
    return: value of Zed according to the pressure value
    '''
    # Parameter evaluation
    Rgas = 0.082057338  # L*atm/K/mol
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    k = (1 + S * (1 - (T / Tcr) ** 0.5)) ** 2
    a = (0.45724 * k * (Rgas * Tcr) ** 2) / Pcr
    b = (0.07780 * Rgas * Tcr) / Pcr
    A = (a * P) / ((Rgas * T) ** 2)
    B = (b * P) / (Rgas * T)
    alpha = - 1 - B + u * B
    beta = A + v * B ** 2 - u * B - u * B ** 2
    gamma = -A * B - v * B ** 2 - v * B ** 3

    # Solver coefficients for Peng-Robinson
    # Objective function: X^3 + p*X + q = 0 --> Zed = X - alpha/3
    p = beta - (alpha ** 2) / 3
    q = 2 * (alpha ** 3) / 27 - alpha * beta / 3 + gamma

    # Discrimant evaluation
    D = (q ** 2) / 4 + (p ** 3) / 27

    # Sign verification:
    if D > 0:
        X = (-q / 2 + (D ** 0.5)) ** (1 / 3) + (-q / 2 - (D ** 0.5)) ** (1 / 3)
        Zed_1 = X - alpha / 3
        Zed_2 = Zed_1
        Zed_3 = Zed_1
    elif D == 0:
        X = (q / 2) ** (1 / 3)
        Zed_1 = -2 * X - alpha / 3
        Zed_2 = X - alpha / 3
        Zed_3 = Zed_2
    elif D < 0:
        r = (-p ** 3 / 27) ** 0.5
        theta = np.arccos(-q / (2 * r))
        Zed_1 = 2 * r ** (1 / 3) * np.cos(theta / 3) - alpha / 3
        Zed_2 = 2 * r ** (1 / 3) * np.cos((2 * mt.pi + theta) / 3) - alpha / 3
        Zed_3 = 2 * r ** (1 / 3) * np.cos((4 * mt.pi + theta) / 3) - alpha / 3

    Zed = [np.real(Zed_1), np.real(Zed_2), np.real(Zed_3)]

    if display == False:
        pass
    elif display == True:
        print('EOS roots:')
        print('Discriminant = ', D)
        print('Zed(1) = ', Zed[0])
        print('Zed(2) = ', Zed[1])
        print('Zed(3) = ', Zed[2])
        print('')

    # Choose the right Zed
    if P / Psat > 1:
        Zed_ = min(Zed)  # Liquid phase
    else:
        Zed_ = max(Zed)  #  Gas phase

    return Zed_


def flash(u, z, Psat, H, P, phi_V_sat, phi_V):
    '''
    sysmple flash system 7x7 for
    thermodynamic state solution
    '''

    # Unknown definition
    y = {'CH4': u[0], 'CO2': u[1], 'H2S': u[2], 'H2O': 1 - u[0] - u[1] - u[2]}  # mol/mol_wet

    x = {'CH4': u[3], 'CO2': u[4], 'H2S': u[5], 'H2O': 1 - u[3] - u[4] - u[5]}  # mol/mol

    alpha = u[6]

    # Equation definition
    eq = np.empty(7)
    eq[0] = P * y['CH4'] * phi_V['CH4'] - phi_V_sat['CH4'] * Psat['CH4'] * x['CH4']
    eq[1] = P * y['CO2'] * phi_V['CO2'] - phi_V_sat['CO2'] * Psat['CO2'] * x['CO2']
    eq[2] = P * y['H2S'] * phi_V['H2S'] - phi_V_sat['H2S'] * Psat['H2S'] * x['H2S']
    eq[3] = P * y['H2O'] * phi_V['H2O'] - phi_V_sat['H2O'] * Psat['H2O'] * x['H2O']
    eq[4] = z['CH4'] - (1 - alpha) * x['CH4'] - alpha * y['CH4']
    eq[5] = z['CO2'] - (1 - alpha) * x['CO2'] - alpha * y['CO2']
    eq[6] = z['H2S'] - (1 - alpha) * x['H2S'] - alpha * y['H2S']

    return eq


def fug_coef_PR(Zed, w, Tcr, Pcr, P, T):
    '''
    evaluation of the fugacity coeffcient
    from the Peng-Robinson EOS in a mixture.
    phi = exp(
      Bi/B*(Zed-1)-A/(2*2^0.5*B)*(Bi/B-2*(Ai/A)^0.5)*
      ln((Zed+B*(1+2^0.5))/(Zed+B*(1-2^0.5)))-ln(Zed-B))
    S = 0.37464 + 1.5422 * w - 0.26992 * w**2
    k = ( 1 + S*(1-(T/Tcr)**0.5) )**2
    a = (0.45724*k*(Rgas*Tcr)**2)/Pcr
    b = (0.07780*Rgas*Tcr)/Pcr
    A = (a*P)/(Rgas*T)**2
    B = (b*P)/(Rgas*T)
    return: value of the fugacity coefficeint accordin to the system pressure
    '''
    # Parameter evaluation
    Rgas = 0.082057338  # L*atm/K/mol
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    k = (1 + S * (1 - (T / Tcr) ** 0.5)) ** 2
    a = (0.45724 * k * (Rgas * Tcr) ** 2) / Pcr
    b = (0.07780 * Rgas * Tcr) / Pcr
    A = (a * P) / ((Rgas * T) ** 2)
    B = (b * P) / (Rgas * T)

    # Parameter gathering
    c1 = Zed - 1
    c2 = A / (2 * (2 ** 0.5) * B)
    c3 = np.log((Zed + B * (1 + 2 ** 0.5)) / (Zed + B * (1 - 2 ** 0.5)))
    c4 = np.log(Zed - B)
    phi = np.exp(c1 - c2 * c3 - c4)

    # check the value
    if np.isnan(phi) == True:
        phi = 0
    else:
        pass
    return phi


def mix_coeff(z, w, Tcr, Pcr, T):
    # Parameter evaluation
    Rgas = 0.082057338  # L*atm/K/mol

    ai = np.empty(shape=(len(z)))
    bi = np.empty(shape=(len(z)))
    S = np.empty(shape=(len(z)))
    k = np.empty(shape=(len(z)))
    zai = np.empty(shape=(len(z)))

    for i in range(len(z)):
        S[i] = 0.37464 + 1.54226 * w[i] - 0.26992 * w[i] ** 2
        k[i] = (1 + S[i] * (1 - (T / Tcr[i]) ** 0.5)) ** 2
        ai[i] = (0.45724 * k[i] * (Rgas * Tcr[i]) ** 2) / Pcr[i]
        bi[i] = (0.07780 * Rgas * Tcr[i]) / Pcr[i]
        zai[i] = z[i] * (ai[i] ** 0.5)

    a = sum(zai) ** 2
    b = sum(bi) / len(bi)

    return [a, b]


def fug_mix_PR(Zed, w, Tcr, Pcr, P, T, amix, bmix):
    '''
    evaluation of the fugacity mixture coeffcient
    from the Peng-Robinson EOS in a mixture.
    phi = exp(
      Bi/B*(Zed-1)-A/(2*2^0.5*B)*(Bi/B-2*(Ai/A)^0.5)*
      ln((Zed+B*(1+2^0.5))/(Zed+B*(1-2^0.5)))-ln(Zed-B))
    Rgas = 0.082057338 # L*atm/K/mol
    S = 0.37464 + 1.5422 * w - 0.26992 * w**2
    k = ( 1 + S*(1-(T/Tcr)**0.5) )**2
    ai = (0.45724*k*(Rgas*Tcr)**2)/Pcr
    bi = (0.07780*Rgas*Tcr)/Pcr
    aij = (ai*aj)*(1-kij) = (sum(zi*(ai)**0.5))**2
    bij = (bi + bj)/2
    A = (a*P)/(Rgas*T)**2
    B = (b*P)/(Rgas*T)
    return: value of the fugacity coefficeint accordin to the system pressure
    '''
    # Parameter evaluation
    Rgas = 0.082057338  # L*atm/K/mol
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    k = (1 + S * (1 - (T / Tcr) ** 0.5)) ** 2
    a = (0.45724 * k * (Rgas * Tcr) ** 2) / Pcr
    b = (0.07780 * Rgas * Tcr) / Pcr
    Amix = (amix * P) / ((Rgas * T) ** 2)
    Bmix = (bmix * P) / (Rgas * T)
    A = (a * P) / ((Rgas * T) ** 2)
    B = (b * P) / (Rgas * T)

    # Parameter gathering
    c1 = B / Bmix * (Zed - 1)
    c2 = A / (2 * (2 ** 0.5) * B)
    c3 = B / Bmix - 2 * (A / Amix) ** 0.5
    c4 = np.log((Zed + Bmix * (1 + 2 ** 0.5)) / (Zed + Bmix * (1 - 2 ** 0.5)))
    c5 = np.log(Zed - Bmix)
    phi_mix = np.exp(c1 + c2 * c3 * c4 - c5)

    # check the value
    if np.isnan(phi_mix) == True:
        phi = 0
    else:
        pass
    return phi_mix


def ln_phi(Z, A, B, phi):
    alpha = -A / (2 * (2 ** 0.5) * B)
    beta_1 = B * (1 + (2 ** 0.5))
    beta_2 = B * (1 - (2 ** 0.5))

    return safe_log(phi) - (Z - 1 - alpha * safe_log((Z + beta_1) / (Z + beta_2)) - safe_log(Z - B))


def fug_mix_helmoltz(specie, T, P, x):
    # importing necessary parameters
    from prop2 import T_c, P_c, w, Ant

    species_vector = list(x.keys())

    z_vec = [x[k] for k in species_vector]
    w_vec = [w[k] for k in species_vector]
    T_c_vec = [T_c[k] for k in species_vector]
    P_c_vec = [P_c[k] / 101325 for k in species_vector]

    x = x[specie]
    w = w[specie]
    T_c = T_c[specie]
    P_c = P_c[specie]

    density_vapour = density(d_vap[specie], T, P, specie)

    phi = VLProperties(6, T, density_vapour, specie)  # ? will be assigned from specie

    P = P / 101325  # Pa to atm for cubic eos
    P_c = P_c / 101325  # Pa to atm for cubic eos
    Psat = Antoine(Ant[specie], T)

    # Parameter evaluation
    Rgas = 0.082057338  # L*atm/K/mol
    S = 0.37464 + 1.54226 * w - 0.26992 * w ** 2
    k = (1 + S * (1 - (T / T_c) ** 0.5)) ** 2
    a = (0.45724 * k * (Rgas * T_c) ** 2) / P_c
    b = (0.07780 * Rgas * T_c) / P_c
    A = (a * P) / ((Rgas * T) ** 2)
    B = (b * P) / (Rgas * T)

    # 1) compressibility factor evaluation from cubic relation through Helmoltz based phi
    Z0 = root_Zed(w, T_c, P_c, P, T, Psat, display=False, u=2, v=-1)  # T[K], P[atm]

    Z = fsolve(ln_phi, Z0, args=(A, B, phi))

    # 2) fugacity mixture coefficient with Z from 1)
    whole_mix = mix_coeff(z_vec, w_vec, T_c_vec, P_c_vec, T)

    mix_whole = {'a': whole_mix[0], 'b': whole_mix[1]}
    phi_mix = fug_mix_PR(Z, w, T_c, P_c, P, T, mix_whole['a'], mix_whole['b'])  # T[K], P[atm]

    phi_cubic = fug_coef_PR(Z0, w, T_c, P_c, P, T)  # T[K], P[atm]
    phi_mix_cubic = fug_mix_PR(Z0, w, T_c, P_c, P, T, mix_whole['a'], mix_whole['b'])  # T[K], P[atm]

    return phi_mix[0]

def dPsatdT(Psat, T):
    """
    Evaluate the rate of change of the vapour pressure
    against the temperature: dPast/dT
    Psat: array
    T: array
    return: rate = incremental ratio of Psat and T
    """
    rate = np.empty(len(Psat))
    for i in range(1, len(Psat)):
        rate[i] = (Psat[i-1] - Psat[i])/(T[i-1] - T[i])
    return rate

def volat_y(y):
    """
    Evaluate the volatility of a species
    in a binary mixture from its composition
    y: array of composition in gas phase
    return: volatility of that species
    """
    return y/(1-y)

def volat_Psat(Psat_A, Psat_B):
    """
    Evaluate the volatility of a species
    in a binary mixture from its vapour pressure
    against the vapour pressure of a referenced one
    Psat_A: vapour pressure of interested species
    Psat_B: vapour pressure of referenced species
    return: volatility of that species
    """
    return Psat_A/Psat_B


# safe log function
def safe_log(x):
    from math import log
    if x <= 0: return   0
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

def mole_day(lks: list, hk: str, Temp_degC, Pres_atm, Height, Diameter, showbool=False):

    dpi = 250

   # evaluate the molar flux of a species ina binary mixture (A-W) where:
    T = Temp_degC  # °C
    p = Pres_atm  # atm
    R = 0.0821e-3 # atm m3 / mol K
    L = Height # m first hypothesis
    d = Diameter # m second hypothesis

    A = lks
    B = hk
    print('molar flux evaluation...')
    print(f'low keys: {A}')
    print(f'high keys: {B}')

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

    # retrieve file path
    # `cwd`: current directory
    dirname = Path.cwd()

    # load resulting compositions (mol/mol)
    file_1 = str(dirname) + '/output/compositions/csv/xy_03atm.csv'
    file_3 = str(dirname) + '/output/compositions/csv/xy_05atm.csv'
    file_5 = str(dirname) + '/output/compositions/csv/xy_10atm.csv'
    file_7 = str(dirname) + '/output/compositions/csv/xy_15atm.csv'



    df_03_atm = pd.read_csv(filepath_or_buffer=file_1, header=0)
    df_05_atm = pd.read_csv(filepath_or_buffer=file_3, header=0)
    df_10_atm = pd.read_csv(filepath_or_buffer=file_5, header=0)
    df_15_atm = pd.read_csv(filepath_or_buffer=file_7, header=0)

    y_H2O = {
        '03atm': df_03_atm['yH2O'],
        '05atm': df_05_atm['yH2O'],
        '10atm': df_10_atm['yH2O'],
        '15atm': df_15_atm['yH2O'],
        }

    # phi = P_H2O(T)/Psat_H2O(T) = p*y_H2O(T,p)/Psat_H2O(T)
    phi = {
        '03atm': [p[0]*y_H2O['03atm'][i]/psat[B][i] for i in range(len(T))],
        '05atm': [p[1]*y_H2O['05atm'][i]/psat[B][i] for i in range(len(T))],
        '10atm': [p[2]*y_H2O['10atm'][i]/psat[B][i] for i in range(len(T))],
        '15atm': [p[3]*y_H2O['15atm'][i]/psat[B][i] for i in range(len(T))],
    }

    # molar fluxes evaluation (mol/m2)
    S = np.pi*d**2/4

    n_03atm ={
        A[0]: [molar_flux(p[0], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['03atm'][i])*S for i in range(len(T))],
        A[1]: [molar_flux(p[0], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['03atm'][i])*S for i in range(len(T))],
        A[2]: [molar_flux(p[0], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['03atm'][i])*S for i in range(len(T))]
    }
    n_05atm ={
        A[0]: [molar_flux(p[1], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['05atm'][i])*S for i in range(len(T))],
        A[1]: [molar_flux(p[1], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['05atm'][i])*S for i in range(len(T))],
        A[2]: [molar_flux(p[1], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['05atm'][i])*S for i in range(len(T))]
    }
    n_10atm ={
        A[0]: [molar_flux(p[2], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['10atm'][i])*S for i in range(len(T))],
        A[1]: [molar_flux(p[2], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['10atm'][i])*S for i in range(len(T))],
        A[2]: [molar_flux(p[2], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['10atm'][i])*S for i in range(len(T))]
    }
    n_15atm ={
        A[0]: [molar_flux(p[3], D[A[0]][i], R, T[i]+273, L, psat[A[0]][i], phi['15atm'][i])*S for i in range(len(T))],
        A[1]: [molar_flux(p[3], D[A[1]][i], R, T[i]+273, L, psat[A[1]][i], phi['15atm'][i])*S for i in range(len(T))],
        A[2]: [molar_flux(p[3], D[A[2]][i], R, T[i]+273, L, psat[A[2]][i], phi['15atm'][i])*S for i in range(len(T))]
    }
    # visualization
    ms = 2
    lw = 0.75

    plt.figure(dpi=dpi)
    filename = 'absolute molar flux [mol/d]'
    plt.plot(T, n_03atm[A[0]], marker='o', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[0]}\,atm$')
    plt.plot(T, n_05atm[A[0]], marker='^', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[1]}\,atm$')
    plt.plot(T, n_10atm[A[0]], marker='s', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[2]}\,atm$')
    plt.plot(T, n_15atm[A[0]], marker='D', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[3]}\,atm$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    plt.ylabel(filename)
    plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.title('CH$_4$')

    plt.figure(dpi=dpi)
    filename = 'absolute molar flux [mol/d]'
    plt.ylabel(filename)
    plt.plot(T, n_03atm[A[1]], marker='o', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[0]}\,atm$')
    plt.plot(T, n_05atm[A[1]], marker='^', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[1]}\,atm$')
    plt.plot(T, n_10atm[A[1]], marker='s', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[2]}\,atm$')
    plt.plot(T, n_15atm[A[1]], marker='D', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[3]}\,atm$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.title('CO$_2$')

    plt.figure(dpi=dpi)
    filename = 'absolute molar flux [mol/d]'
    plt.ylabel(filename)
    plt.plot(T, n_03atm[A[2]],marker='o', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[0]}\,atm$')
    plt.plot(T, n_05atm[A[2]],marker='^', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[1]}\,atm$')
    plt.plot(T, n_10atm[A[2]],marker='s', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[2]}\,atm$')
    plt.plot(T, n_15atm[A[2]],marker='D', markersize=ms, linewidth=lw, label=f'P: ${Pres_atm[3]}\,atm$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    plt.title('H$_2$S')
    plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    plt.figure(dpi=dpi)
    filename = 'absolute molar flux [mol/d]'
    plt.ylabel(filename)
    plt.plot(T, n_03atm[A[0]]/np.max(n_03atm[A[0]]),marker='o', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[0]}\,atm$')
    plt.plot(T, n_05atm[A[0]]/np.max(n_05atm[A[0]]),marker='^', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[1]}\,atm$')
    plt.plot(T, n_10atm[A[0]]/np.max(n_10atm[A[0]]),marker='s', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[2]}\,atm$')
    plt.plot(T, n_15atm[A[0]]/np.max(n_15atm[A[0]]),marker='D', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[3]}\,atm$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.title('CH$_4$')

    plt.figure(dpi=dpi)
    filename = 'absolute molar flux [mol/d]'
    plt.ylabel(filename)
    plt.plot(T, n_03atm[A[1]]/np.max(n_03atm[A[1]]),marker='o', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[0]}\,atm$')
    plt.plot(T, n_05atm[A[1]]/np.max(n_05atm[A[1]]),marker='^', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[1]}\,atm$')
    plt.plot(T, n_10atm[A[1]]/np.max(n_10atm[A[1]]),marker='s', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[2]}\,atm$')
    plt.plot(T, n_15atm[A[1]]/np.max(n_15atm[A[1]]),marker='D', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[3]}\,atm$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    plt.xlabel('temperature [°C]');
    plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.title('CO$_2$')

    plt.figure(dpi=dpi)
    filename = 'absolute molar flux [mol/d]'
    plt.ylabel(filename)
    plt.plot(T, n_03atm[A[2]]/np.max(n_03atm[A[2]]),marker='o', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[0]}\,atm$')
    plt.plot(T, n_05atm[A[2]]/np.max(n_05atm[A[2]]),marker='^', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[1]}\,atm$')
    plt.plot(T, n_10atm[A[2]]/np.max(n_10atm[A[2]]),marker='s', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[2]}\,atm$')
    plt.plot(T, n_15atm[A[2]]/np.max(n_15atm[A[2]]),marker='D', linewidth=lw, markersize=ms, label=f'P: ${Pres_atm[3]}\,atm$')
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    plt.xlabel('temperature [°C]');
    plt.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    plt.title('H$_2$S')

    if showbool == True: plt.show()
    else: pass

    return [n_03atm, n_05atm, n_10atm, n_15atm]

