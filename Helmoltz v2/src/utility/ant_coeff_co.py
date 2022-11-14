import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

coeff = np.array([6.24020, 230.270, 260.010])


# functions
def mmhg_to_bar(x): return x * 0.00133322


def degc_to_k(x): return x + 275.15


def antonie_expr(x, A, B, C): return 10 ** (A - B / (x + C))


def lin_ant(x, y, A, B, C): return x*(1/A-1/np.log10(y))*A/(A*C-B) + C/(A*C-B)


# data preparation
p_mmhg = np.empty(50)
inv_logp_mmhg = np.empty(50)
inv_logp_bar = np.empty(50)
p_bar = np.empty(50)
T_k = np.empty(50)
T_degc = np.linspace(-210, -165, 50)

for i in range(len(T_degc)):
    p_mmhg[i] = antonie_expr(T_degc[i], coeff[0], coeff[1], coeff[2])
    p_bar[i] = mmhg_to_bar(p_mmhg[i])
    T_k[i] = degc_to_k(T_degc[i])
    inv_logp_mmhg[i] = 1/np.log10(p_mmhg[i])
    inv_logp_bar[i] = 1/np.log10(p_bar[i])


# fitting
plt.figure()
plt.plot(T_degc, lin_ant(T_degc, p_mmhg, *coeff), 'ko')
plt.xlabel('T_degC')
plt.ylabel('1/log10(P)_mmhg')
# plt.show()
