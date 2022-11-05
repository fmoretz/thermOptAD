from prop import *
from func import *
import math as mt
from scipy.optimize import fsolve
    
def ln_phi(Z, A, B, phi):
    
    alpha  = -A/(2*(2**0.5)*B)
    beta_1 = B*(1+(2**0.5)) 
    beta_2 = B*(1-(2**0.5)) 
    
    return mt.log(phi) - ( Z - 1 - alpha * mt.log( (Z+beta_1)/(Z+beta_2) ) - mt.log(Z-B) ) 

def fug_mix_helmoltz(specie):
    
    # importing necessary parameters
    from prop import T_c, P_c, w, Ant, R_m
    from initial import T, P, z
    
    species_vector = list(z.keys())
    
    z_vec   = [z[k]          for k in species_vector]
    w_vec   = [w[k]          for k in species_vector]
    T_c_vec = [T_c[k]        for k in species_vector]
    P_c_vec = [P_c[k]/101325 for k in species_vector]
    
    z   = z[specie]
    w   = w[specie]
    T_c = T_c[specie]
    P_c = P_c[specie]
    
    density_vapour = density(d_vap[specie], T, P, specie)

    phi = VLProperties(6, T, density_vapour, specie)

    P   = P   / 101325             # Pa to atm for cubic eos relations
    P_c = P_c / 101325             # Pa to atm for cubic eos relations
    Psat = Antoine(Ant[specie], T) 

    # Parameter evaluation
    Rgas = 0.082057338 # L*atm/K/mol
    S = 0.37464 + 1.54226 * w - 0.26992 * w**2
    k = ( 1 + S*(1-(T/T_c)**0.5) )**2
    a = (0.45724*k*(Rgas*T_c)**2)/P_c
    b = (0.07780*Rgas*T_c)/P_c
    A = (a*P)/((Rgas*T)**2)
    B = (b*P)/(Rgas*T)

    # 1) compressibility factor evaluation from cubic relation through Helmoltz based phi
    Z0 = root_Zed(w, T_c, P_c, P, T, Psat, display = False, u = 2, v = -1) # T[K], P[atm]

    rho = density(d_vap[specie], T, P, specie)
    
    d_r = rho/r_c[specie]
    tau = T_c[specie]/T
    
    Z_rho = dphi_r(d_r,tau,species) # P / (R_m[specie] * T * rho)
    Z = fsolve(ln_phi, Z0, args = (A, B, phi))

    # 2) fugacity mixture coefficient with Z from 1)
    whole_mix = mix_coeff(z_vec, w_vec, T_c_vec, P_c_vec, T)

    mix_whole = {
        'a': whole_mix[0],
        'b': whole_mix[1]
    }
    phi_mix = fug_mix_PR(Z, w, T_c, P_c, P, T, mix_whole['a'], mix_whole['b']) # T[K], P[atm]

    phi_cubic     = fug_coef_PR(Z0, w, T_c, P_c, P, T) # T[K], P[atm]
    phi_mix_cubic = fug_mix_PR(Z0, w, T_c, P_c, P, T, mix_whole['a'], mix_whole['b']) # T[K], P[atm]

    print('# ======================================= #')
    print('Comparison:')
    print(f'Z from Helmoltz: {Z[0]}')
    print(f'Z from Cubic PR: {Z0}')
    print(f'Z from density : {Z_rho}')
    print('')
    print(f'fug.pure from Helmoltz: {phi}')
    print(f'fug.pure from Cubic PR: {phi_cubic}')
    print('')
    print(f'fug.mix from Helmoltz: {phi_mix[0]}')
    print(f'fug.mix from Cubic PR: {phi_mix_cubic}')
    print('# ======================================= #')
    
    return phi_mix[0]

specie = 'H2O'
fug_mix_helmoltz(specie)