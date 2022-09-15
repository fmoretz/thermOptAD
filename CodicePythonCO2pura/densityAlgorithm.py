#---------------------------------------------------------
# This module contains the algorithms developed in order
# to obtain carbon dioxide density which is necessary for
# the calculation of phase equilibrium conditions and of
# properties from Span & Wagner's Helmholtz energy model
#---------------------------------------------------------

from scipy import optimize as opt
# Module containing the derivatives of the Helmholtz energy
import SWderivativeFunctions as SWdF

def density(m,T,P):
    """
    #----------------------------------------------
    # Function providing density of carbon dioxide
    # at a given pressure and temperature
    #----------------------------------------------
    """

    # Set up start value of density in units of Kg*m^(-3)
    if m==1:
        x0 = 1200.0
    elif m==2:
        x0 = 4.0

    #---------------------------------------------
    # Newton-Raphson algorithm used to calculate
    # the zero of a non-analytic function
    #---------------------------------------------

    x = opt.newton(rho,x0,fprime=drho, args=(T,P), tol = 1.0E-6, maxiter=100)

    return x

def rho(x,T,P):
    """
    #---------------------------------------------------
    # Function describing relation between pressure and
    # density: P = RT*rho*(1 + rho_ridotta*(phi_r)')
    #---------------------------------------------------
    """

    # Set up reduced density and the inverse reduced temperature
    d_r  = x/r_c
    tau  = T_c/T

    # Derivative of the residual part of the Helmholtz energy with respect to reduced density
    dphi = SWdF.dphi_r(d_r,tau)

    #-----------------------------------------------------
    # Function rho(x) --> variable x = density [Kg/m^(3)]
    #-----------------------------------------------------
    return x*R_m*T*(1 + d_r*dphi) - P

def drho(x,T,P):
    """
    #-----------------------------------------------
    # Function describing the derivative of rho(x)
    #-----------------------------------------------
    """

    # Set up reduced density and the inverse reduced temperature
    d_r  = x/r_c
    tau  = T_c/T

    # Derivatives of the residual part of the Helmholtz energy with respect to reduced density
    dphi    = SWdF.dphi_r(d_r,tau)
    ddphi_r = SWdF.d2phi_r(d_r,tau)

    #------------------------------------------------------
    # Function rho'(x) --> variable x = density [Kg/m^(3)]
    #------------------------------------------------------
    return R_m*T*(1 + 2*d_r*dphi + (d_r**2)*ddphi_r)
