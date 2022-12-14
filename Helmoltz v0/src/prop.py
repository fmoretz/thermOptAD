import numpy as np
# Parameters describing the critical point of the compunds

C = {# A, B, C, D, E, and F
    'CH4': [3.135E+1, -1.308E+3, 0, -3.261E+0, 2.942E-5, 2],
    'CO2': [1.336E+2, -4.735E+3, 0, -2.127E+1, 4.091E-2, 1],
    'H2S': [7.868E+1, -3.840E+3, 0, -1.120E+1, 1.885E-2, 1],
    'H2O': [6.593E+1, -7.228E+3, 0, -7.177E+0, 4.031E-6, 2]
}

Ant = {# A, B, C
    'CH4': [3.98950,    443.028,	 -0.49],
    'CO2': [6.81228,	1301.679,	-3.494],
    'H2S': [4.52887,	958.5870,	-0.539],
    'H2O': [6.20963,	2354.731,	 7.559]
}


kH = {# A, B, C and D
    'CH4': [-142.234, 26.5960, 32.308, -0.090],
    'CO2': [-183.691, 676.278, 39.068, -0.098],
    'H2S': [-186.884, 1090.12, 39.011, -0.094]
}

w = {
    'CH4': 0.011,
    'CO2': 0.239,
    'H2S': 0.081,
    'H2O': 0.344
} # --

T_c = {
    'CH4': 190.564,
    'CO2': 304.1282,
    'H2S': 373.1,
    'H2O': 647.096
} # K

P_c = {
    'CH4': 4592200,
    'CO2': 7377300,
    'H2S': 9000000,
    'H2O': 22064000
} # Pa

r_c = {
    'CH4': 126.66,
    'CO2': 467.60,
    'H2S': 347.63,
    'H2O': 322
} # Kg*m^(-3)

# Specific gas constant
R_m = {
    'CH4': 518.267548,
    'CO2': 188.924269,
    'H2S': 243.961815,
    'H2O': 461.51805
} # J*kg^(-1)*K^(-1)
                 
# Universal gas constant
R = 8.31446262               # In units of J*mol^(-1)*K^(-1)

# Start values for the liquid density
d_liq = {
    'CH4': 500,
    'CO2': 1200,
    'H2S': 1000,
    'H2O': 1200
} # Kg*m^(-3)

# Start values for the vapor density
d_vap = {
    'CH4': 0.01,
    'CO2': 4.0,
    'H2S': 1.0,
    'H2O': 1200 # 0.01
} # Kg*m^(-3)

# Parameters of the ideal-gas part of the Helmholtz energy
a = {# a_1, a_2, a_3 and a forth term
    'CH4': [9.91243972, -6.33270087, 3.0016, 0],
    'CO2': [8.37304456, -3.70454304, 2.5000, 0],
    'H2S': [-4.0740770957, 3.7632137341, 3, (0.14327e-5)/3.75],
    'H2O': [-8.32044648201, 6.6832105268, 3.00632, 0]
}

a0 = {# a_4, a_5, a_6, a_7 and a_8
    'CH4': [0.008449, 4.6942, 3.4865, 1.6572, 1.4115],
    'CO2': [1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678],
    'H2S': [1.1364, 1.9721, 1e-10, 1e-10],
    'H2O': [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
}

th0 = {# theta_4, theta_5, theta_6, theta_7 and theta_8
    'CH4': [3.4004324, 10.26951575, 20.43932747, 29.93744884, 79.13351945],
    'CO2': [3.15163, 6.11190, 6.77708, 11.32384, 27.08792],
    'H2S': [1823.0/373.1, 3965.0/373.1, 1e-10, 1e-10],
    'H2O': [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]
}

# Parameters of the residual part of the Helmholtz energy
#---------------------------First term--------------------------------------------------------------------------------------------------------
n1 = {
    'CH4': [0.04367901028, 0.6709236199, -1.765577859, 0.8582330241, -1.206513052, 0.512046722, -0.04000010791e-2, -0.01247842423, 
            0.03100269701, 0.1754748522e-2, -0.3171921605e-5, -0.224034684e-5, 0.2947056156e-6],
    'CO2': [0.38856823203161,2.938547594274,-5.5867188534934,-0.76753199592477,0.31729005580416,0.54803315897767,0.12279411220335],
    'H2S': [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.00019076],
    'H2O': [0.012533547935523, 7.8957634722828, -8.7803203303561, 0.31802509345418, -0.26145533859358, -0.78199751687981e-2, 0.88089493102134e-2]
}

d1 = {
    'CH4': [1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 8, 9, 10],
    'CO2': [1, 1, 1, 1, 2, 2, 3],
    'H2S': [1, 1, 1, 2, 3, 7],
    'H2O': [1, 1, 1, 2, 2, 3, 4]
}

t1 = {
    'CH4': [-0.5, 0.5, 1, 0.5, 1, 1.5, 4.5, 0, 1, 3, 1, 3, 3],
    'CO2': [0.00, 0.75, 1.00, 2.00, 0.75, 2.00, 0.75],
    'H2S': [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],
    'H2O': [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1.0]
}
#---------------------------Second term------------------------------------------------------------------------------------------------------------------
n2 = {
    'CH4': [0.1839487909, 0.1511883679, -0.4289363877, 0.06894002446, -0.01408313996, -0.0306305483, -0.02969906708, -0.01932040831, -0.1105739959, 
            0.09952548995, 0.8548437825e-2, -0.06150555662, -0.04291792423, -0.0181320729, 0.0344590476, -0.238591945e-2, -0.01159094939, 0.06641693602,
            -0.0237154959, -0.03961624905, -0.01387292044, 0.03389489599, -0.2927378753e-2],
    'CO2': [2.1658961543220, 1.5841735109724,-0.23132705405503,5.8116916431436e-2,-0.55369137205382, 0.48946615909422,-2.4275739843501e-2,
            0.062494790501678,-0.12175860225246, -0.37055685270086,-1.6775879700426e-2,-0.11960736637987,-0.045619362508778, 3.5612789270346e-2,
            -7.4427727132052e-3,-1.7395704902432e-3,-0.021810121289527,0.024332166559236,-3.7440133423463e-2,0.14338715756878,-0.13491969083286,
            -0.02315122505348, 1.2363125492901e-2,0.002105832197294,-3.3958519026368e-4, 5.5993651771592e-3, -3.0335118055646e-4],
    'H2S': [0.20227, -0.0045348, -0.22230, -0.034714, -0.014885, 0.0074154],
    'H2O': [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002, -0.25709043003438, 0.16074868486251, -0.040092828925807,
            0.39343422603254e-6, -0.75941377088144e-5, 0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8, 0.36582165144204e-6, 
            -0.13251180074668e-11, -0.62639586912454e-9, -0.10793600908932, 0.017611491008752, 0.22132295167546, -0.40247669763528, 0.58083399985759,
            0.49969146990806e-2, -0.031358700712549, -0.74315929710341, 0.4780732991548, 0.020527940895948, -0.13636435110343, 0.014180634400617,
            0.83326504880713e-2, -0.029052336009585, 0.038615085574206, -0.020393486513704, -0.16554050063734e-2, 0.19955571979541e-2, 0.15870308324157e-3,
            0.1638856834253e-4, 0.043613615723811, 0.034994005463765, -0.076788197844621, 0.022446277332006, -0.62689710414685e-4, -0.55711118565645e-9,
            -0.19905718354408, 0.31777497330738, -0.11841182425981]
}

d2 = {
    'CH4': [1, 1, 1, 2, 4, 5, 6, 1, 2, 3, 4, 4, 3, 5, 5, 8, 2, 3, 4, 4, 4, 5, 6],
    'CO2': [1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8],
    'H2S': [2, 5, 1, 4, 3, 4],
    'H2O': [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
}

t2 = {
    'CH4': [0, 1, 2, 0, 0, 2, 2, 5, 5, 5, 2, 4, 12, 8, 10, 10, 10, 14, 12, 18, 22, 18, 14],
    'CO2': [1.5,1.5,2.5,0.0,1.5,2.0,0.0,1.0,2.0,3.0,6.0,3.0,6.0,8.0,6.0,0.0,7.0,12.0,16.0,22.0,24.0,16.0,24.0,8.0,2.0,28.0,14.0],
    'H2S': [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
    'H2O': [4,6,12,1,5,4,2,13,9,3,4,11,4,13,1,7,1,9,10,10,3,7,10,10,6,10,10,1,2,3,4,8,6,9,8,16,22,23,23,10,50,44,46,50]
}

c2 = {
    'CH4': [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4],
    'CO2': [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6],
    'H2S': [1, 1, 2, 2, 3, 3],
    'H2O': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]
}
#---------------------------Third term---------------------------------------------------------------------------------------------------------
n3 = {
    'CH4': [0.9324799946e-4, -6.287171518, 12.71069467, -6.423953466],
    'CO2': [-213.6548868832,26641.569149272,-24027.212204557,-283.41603423999,212.47284400179],
    'H2S': [0, 0, 0, 0],
    'H2O': [-31.306260323435, 31.546140237781, -2521.3154341695]
}

d3 = {
    'CH4': [2, 0, 0, 0],
    'CO2': [2, 2, 2, 3, 3],
    'H2S': [0, 0, 0, 0],
    'H2O': [3, 3, 3]
}

t3 = {
    'CH4': [2, 0, 1, 2],
    'CO2': [1.0, 0.0, 1.0, 3.0, 3.0],
    'H2S': [0, 0, 0, 0],
    'H2O': [0, 1, 4]
}

alpha3 = {
    'CH4': [20, 40, 40, 40],
    'CO2': [25, 25, 25, 15, 20],
    'H2S': [0, 0, 0, 0],
    'H2O': [20, 20, 20]
}

beta3 = {
    'CH4': [200, 250, 250, 250],
    'CO2': [325, 300, 300, 275, 275],
    'H2S': [0, 0, 0, 0],
    'H2O': [150, 150, 250]
}

gamma3 = {
    'CH4': [1.07, 1.11, 1.11, 1.11],
    'CO2': [1.16, 1.19, 1.19, 1.25, 1.22],
    'H2S': [0, 0, 0, 0],
    'H2O': [1.21, 1.21, 1.25]
}
#---------------------------Fourth term---------------------------------------------------------------------------------------------------------
n4 = {
    'CH4': [0, 0, 0],
    'CO2': [-0.66642276540751, 0.72608632349897, 0.055068668612842],
    'H2S': [0, 0, 0],
    'H2O': [-0.14874640856724, 0.31806110878444]
}

C4 = {
    'CH4': [0, 0, 0],
    'CO2': [10.0, 10.0, 12.5],
    'H2S': [0, 0, 0],
    'H2O': [28, 32]
}

alpha4 = {
    'CH4': [0, 0, 0],
    'CO2': [3.5, 3.5, 3.0],
    'H2S': [0, 0, 0],
    'H2O': [3.5, 3.5]
}

beta4 = {
    'CH4': [0, 0, 0],
    'CO2': [0.3, 0.3, 1.0],
    'H2S': [0, 0, 0],
    'H2O': [0.2, 0.2]
}

bi4 = {
    'CH4': [0, 0, 0],
    'CO2': [0.875, 0.875, 0.875],
    'H2S': [0, 0, 0],
    'H2O': [0.85, 0.95]
}

Di4 = {
    'CH4': [0, 0, 0],
    'CO2': [275, 275, 275],
    'H2S': [0, 0, 0],
    'H2O': [700, 800]
}

Ai4 = {
    'CH4': [0, 0, 0],
    'CO2': [0.7, 0.7, 0.7],
    'H2S': [0, 0, 0],
    'H2O': [0.32, 0.32]
}