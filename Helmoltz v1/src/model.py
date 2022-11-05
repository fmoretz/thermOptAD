import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sympy import sympify
from func2 import *
import os
from initial import *


def model(Temperature, Pressure, species):
    plt.figure(num=None, dpi=80, facecolor='w', edgecolor='k')
    matplotlib.rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    formatter.set_useOffset(0)

    alpha0 = z['CH4'] + z['CO2'] + z['H2S']

    alpha_TP = np.empty(shape=(len(Pressure), len(Temperature)))
    x_CH4 = np.empty(shape=(len(Pressure), len(Temperature)))
    y_CH4 = np.empty(shape=(len(Pressure), len(Temperature)))
    x_CO2 = np.empty(shape=(len(Pressure), len(Temperature)))
    y_CO2 = np.empty(shape=(len(Pressure), len(Temperature)))
    x_H2S = np.empty(shape=(len(Pressure), len(Temperature)))
    y_H2S = np.empty(shape=(len(Pressure), len(Temperature)))
    x_H2O = np.empty(shape=(len(Pressure), len(Temperature)))
    y_H2O = np.empty(shape=(len(Pressure), len(Temperature)))
    Vap_TP = np.empty(shape=(len(Pressure), len(Temperature)))
    Liq_TP = np.empty(shape=(len(Pressure), len(Temperature)))

    i = 0  # index for cycling in T for-loop
    j = 0  # index for cycling in P for-loop

    for T in Temperature:
        T = T + 273  # K

        for P in Pressure:
            P = P * 101325

            # initialization
            Ar_V = {key: 0 for key in species}  # --
            entropy_V = {key: 0 for key in species}  # J/mol/K
            enthalpy_V = {key: 0 for key in species}  # J/mol/K
            cp_V = {key: 0 for key in species}  # J/mol/K
            cv_V = {key: 0 for key in species}  # J/mol/K
            cubicExp_V = {key: 0 for key in species}  # --
            fuga_V = {key: 0 for key in species}  # --
            fuga_m = {key: 0 for key in species}  # --
            Ar_L = {key: 0 for key in species}  # --
            entropy_L = {key: 0 for key in species}  # J/mol/K
            enthalpy_L = {key: 0 for key in species}  # J/mol/K
            cp_L = {key: 0 for key in species}  # J/mol/K
            cv_L = {key: 0 for key in species}  # J/mol/K
            cubicExp_L = {key: 0 for key in species}  # --
            fuga_L = {key: 0 for key in species}  # --
            Psat = {key: 0 for key in species}  # Pa
            dens_L = {key: 0 for key in species}  # kg * m^(-3)
            dens_V = {key: 0 for key in species}  # kg * m^(-3)
            x = {key: 0 for key in species}  # mol/mol
            y = {key: 0 for key in species}  # mol/mol
            Lx = {key: 0 for key in species}  # kmol/d
            Vy = {key: 0 for key in species}  # kmol/d

            for specie in species:

                if specie == 'CO':
                    Psat[specie] = Antoine(Ant[specie], T-273.15) * 133.322  # mmHg to Pa
                    
                if specie == 'CH4':
                    Psat[specie] = Henry(kH[specie], T) * 1e3 # kPa to Pa
                    
                if specie == 'CO2':
                    Psat[specie] = Henry(kH[specie], T) * 1e3 # kPa to Pa
                    
                else:
                    Psat[specie] = Antoine(Ant[specie], T) * 101325 # kPa to Pa
                
                dens_L[specie] = density(d_liq[specie], T, Psat[specie], specie)
                dens_V[specie] = density(d_vap[specie], T, P, specie)
                Ar_L[specie] = sympify(VLProperties(0, T, dens_L[specie], specie))
                entropy_L[specie] = sympify(VLProperties(1, T, dens_L[specie], specie))
                enthalpy_L[specie] = sympify(VLProperties(2, T, dens_L[specie], specie))
                cp_L[specie] = sympify(VLProperties(3, T, dens_L[specie], specie))
                cv_L[specie] = sympify(VLProperties(4, T, dens_L[specie], specie))
                cubicExp_L[specie] = sympify(VLProperties(5, T, dens_L[specie], specie))
                fuga_L[specie] = sympify(VLProperties(6, T, dens_L[specie], specie))
                fuga_m[specie] = fug_mix_helmoltz(specie, T, P)
                Ar_V[specie] = sympify(VLProperties(0, T, dens_V[specie], specie))
                entropy_V[specie] = sympify(VLProperties(1, T, dens_V[specie], specie))
                enthalpy_V[specie] = sympify(VLProperties(2, T, dens_V[specie], specie))
                cp_V[specie] = sympify(VLProperties(3, T, dens_V[specie], specie))
                cv_V[specie] = sympify(VLProperties(4, T, dens_V[specie], specie))
                cubicExp_V[specie] = sympify(VLProperties(5, T, dens_V[specie], specie))
                fuga_V[specie] = sympify(VLProperties(6, T, dens_V[specie], specie))

            k = {
                key: Psat[key] * fuga_L[key] / (P * fuga_m[key]) for key in species
            } # Pa/Pa

            # Rachford-Rice Solution
            reduced_z = {key: z[key] for key in species}
            reduced_k = {key: k[key] for key in species}
            
            alpha = fsolve(
                func=RR,
                x0=alpha0,
                args=(
                    list(reduced_z.values()),
                    list(reduced_k.values())
                )
            )

            alpha = alpha[0]
            V = alpha * F
            L = F - V

            for specie in species:
                x[specie] = z[specie] / (1 + alpha * (k[specie] - 1))
                y[specie] = x[specie] * k[specie]
                Vy[specie] = V * y[specie]
                Lx[specie] = L * x[specie]

            # array allocation of pressure(j)-temperature(i) dependent variables
            alpha_TP[j][i] = alpha
            x_CH4[j][i] = x['CH4']
            y_CH4[j][i] = y['CH4']
            x_CO2[j][i] = x['CO2']
            y_CO2[j][i] = y['CO2']
            x_H2S[j][i] = x['H2S']
            y_H2S[j][i] = y['H2S']
            x_H2O[j][i] = x['H2O']
            y_H2O[j][i] = y['H2O']
            Vap_TP[j][i] = V
            Liq_TP[j][i] = L

            # Saving thermodynamic properties in a txt file
            P_str = str(P / 101325)
            P_str = P_str.replace('.', '')

            filename = f'results_{int(T - 273)}C_' + P_str + 'atm'
            # cd = os.chdir('')
            path = '/Users/fmoretta/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Thermodynamic Optimization ' \
                   'AD/Helmoltz v0/output/thermo parameters/'

            file = path + filename + '.txt'
            
            with open(file, 'w') as res:
                res.write('# ======= Thermodynamic parameters for species and conditions ======= #')
                res.write('\nOperative Conditions:')
                res.write(f'\nTemperature____: {"{:.3f}".format(T - 273)} [°C]')
                res.write(f'\nPressure_______: {"{:.3f}".format(P / 101325)} [atm]')
                res.write(f'\nVapour fraction: {"{:.3f}".format(alpha)} [--]')
                res.write('\n')
                for specie in species:
                    res.write(f'\nSpecies: {specie} -- phase: Liquid')
                    res.write(f"\nDimensionless SW's Helmholtz energy: {'{:.5f}'.format(Ar_L[specie])} [--]")
                    res.write(f'\nMolar entropy______________________: {"{:.2f}".format(entropy_L[specie])} [J/mol/K]')
                    res.write(f'\nMolar enthalpy_____________________: {"{:.2f}".format(enthalpy_L[specie])} [J/mol/K]')
                    res.write(f'\nIsobaric heat capacity_____________: {"{:.3f}".format(cp_L[specie])} [J/mol/K]')
                    res.write(f'\nIsochoric heat capacity____________: {"{:.3f}".format(cv_L[specie])} [J/mol/K]')
                    res.write(f'\nCubic expansion coefficient________: {"{:.5f}".format(cubicExp_L[specie])} [1/K]')
                    res.write(f'\nFugacity coeffcient________________: {"{:.5}".format(fuga_L[specie])} [--]')
                    res.write('\n')
                    res.write(f'\nSpecies: {specie} -- phase: Vapour')
                    res.write(f"\nDimensionless SW's Helmholtz energy: {'{:.5f}'.format(Ar_V[specie])} [--]")
                    res.write(f'\nMolar entropy______________________: {"{:.2f}".format(entropy_V[specie])} [J/mol/K]')
                    res.write(f'\nMolar enthalpy_____________________: {"{:.2f}".format(enthalpy_V[specie])} [J/mol/K]')
                    res.write(f'\nIsobaric heat capacity_____________: {"{:.3f}".format(cp_V[specie])} [J/mol/K]')
                    res.write(f'\nIsochoric heat capacity____________: {"{:.3f}".format(cv_V[specie])} [J/mol/K]')
                    res.write(f'\nCubic expansion coefficient________: {"{:.5f}".format(cubicExp_V[specie])} [1/K]')
                    res.write(f'\nFugacity coeffcient________________: {"{:.5f}".format(fuga_V[specie])} [--]')
                    res.write(f'\nFugacity coeffcient (mixture)______: {"{:.5f}".format(fuga_m[specie])} [--]')
                    res.write('\n')
                    res.write('\n')
                res.close()

            if j == 7:
                j = 0
            else:
                j = j + 1
        i = i + 1

    labels = [
        'P = ' + str(Pressure[0]) + ' atm',
        'P = ' + str(Pressure[1]) + ' atm',
        'P = ' + str(Pressure[2]) + ' atm',
        'P = ' + str(Pressure[3]) + ' atm',
        'P = ' + str(Pressure[4]) + ' atm',
        'P = ' + str(Pressure[5]) + ' atm',
        'P = ' + str(Pressure[6]) + ' atm',
        'P = ' + str(Pressure[7]) + ' atm'
    ]

    # alpha plot
    fig, axes = plt.subplots(1, 1, figsize=(10, 10))

    axes.plot(Temperature, alpha_TP[0], 'k', linewidth=2, label=labels[0])
    axes.plot(Temperature, alpha_TP[1], 'b', linewidth=2, label=labels[1])
    axes.plot(Temperature, alpha_TP[2], 'c', linewidth=2, label=labels[2])
    axes.plot(Temperature, alpha_TP[3], 'm', linewidth=2, label=labels[3])
    axes.plot(Temperature, alpha_TP[4], 'r', linewidth=2, label=labels[4])
    axes.plot(Temperature, alpha_TP[5], 'orange', linewidth=2, label=labels[5])
    axes.plot(Temperature, alpha_TP[6], 'brown', linewidth=2, label=labels[6])
    axes.plot(Temperature, alpha_TP[7], 'lime', linewidth=2, label=labels[7])
    axes.yaxis.set_major_formatter(formatter)
    axes.set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes.set_ylabel(r"$\alpha\,[\mathrm{--}]$", fontsize=15)

    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes.set_xlim(lim_min, lim_max)
    axes.xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    lim_min = min(alpha_TP[7]) * 0 + 0.008
    lim_max = 0.02
    axes.set_ylim(lim_min, lim_max)
    axes.yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    axes.grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    axes.legend(loc='upper left', fontsize=12)

    path_for_figures = '/Users/fmoretta/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Thermodynamic Optimization ' \
                       'AD/Helmoltz v1/output//VLE figures/'
    filename = path_for_figures + 'alpha_TP_plot.svg'
    fig.savefig(filename, dpi=800)

    # V-L plot
    fig2, axes2 = plt.subplots(1, 2, figsize=(20, 10))

    axes2[0].plot(Temperature, Vap_TP[0], 'k', linewidth=2, label=labels[0])
    axes2[0].plot(Temperature, Vap_TP[1], 'b', linewidth=2, label=labels[1])
    axes2[0].plot(Temperature, Vap_TP[2], 'c', linewidth=2, label=labels[2])
    axes2[0].plot(Temperature, Vap_TP[3], 'm', linewidth=2, label=labels[3])
    axes2[0].plot(Temperature, Vap_TP[4], 'r', linewidth=2, label=labels[4])
    axes2[0].plot(Temperature, Vap_TP[5], 'orange', linewidth=2, label=labels[5])
    axes2[0].plot(Temperature, Vap_TP[6], 'brown', linewidth=2, label=labels[6])
    axes2[0].plot(Temperature, Vap_TP[7], 'lime', linewidth=2, label=labels[7])

    axes2[0].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes2[0].set_ylabel(r"$Vapour$ $flow\,[\mathrm{kmol/d}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature);
    lim_max = max(Temperature)
    axes2[0].set_xlim(lim_min, lim_max)
    axes2[0].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(Vap_TP[7]) / 1.02;
    lim_max = max(Vap_TP[0])
    axes2[0].set_ylim(lim_min, lim_max)
    axes2[0].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    axes2[0].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    axes2[0].legend(loc='upper left', fontsize=12)

    axes2[1].plot(Temperature, Liq_TP[0], 'k', linewidth=2, label=labels[0])
    axes2[1].plot(Temperature, Liq_TP[1], 'b', linewidth=2, label=labels[1])
    axes2[1].plot(Temperature, Liq_TP[2], 'c', linewidth=2, label=labels[2])
    axes2[1].plot(Temperature, Liq_TP[3], 'm', linewidth=2, label=labels[3])
    axes2[1].plot(Temperature, Liq_TP[4], 'r', linewidth=2, label=labels[4])
    axes2[1].plot(Temperature, Liq_TP[5], 'orange', linewidth=2, label=labels[5])
    axes2[1].plot(Temperature, Liq_TP[6], 'brown', linewidth=2, label=labels[6])
    axes2[1].plot(Temperature, Liq_TP[7], 'lime', linewidth=2, label=labels[7])

    axes2[1].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes2[1].set_ylabel(r"$Liquid$ $flow\,[\mathrm{kmol/d}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature);
    lim_max = max(Temperature)
    axes2[1].set_xlim(lim_min, lim_max)
    axes2[1].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(Liq_TP[0]) - 5;
    lim_max = max(Liq_TP[7]) + 5
    axes2[1].set_ylim(lim_min, lim_max)
    axes2[1].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    axes2[1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    axes2[1].legend(loc='upper left', fontsize=12)

    path_for_figures = '/Users/fmoretta/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Thermodynamic Optimization ' \
                       'AD/Helmoltz v1/output/VLE figures/'
    filename = path_for_figures + 'V_L_TP_plot.svg'
    fig2.savefig(filename, dpi=800)

    # x-y plot
    fig3, axes3 = plt.subplots(4, 2, figsize=(15, 20))

    axes3[0, 0].plot(Temperature, x_CH4[0], 'k', linewidth=2, label=labels[0])
    axes3[0, 0].plot(Temperature, x_CH4[1], 'b', linewidth=2, label=labels[1])
    axes3[0, 0].plot(Temperature, x_CH4[2], 'c', linewidth=2, label=labels[2])
    axes3[0, 0].plot(Temperature, x_CH4[3], 'm', linewidth=2, label=labels[3])
    axes3[0, 0].plot(Temperature, x_CH4[4], 'r', linewidth=2, label=labels[4])
    axes3[0, 0].plot(Temperature, x_CH4[5], 'orange', linewidth=2, label=labels[5])
    axes3[0, 0].plot(Temperature, x_CH4[6], 'brown', linewidth=2, label=labels[6])
    axes3[0, 0].plot(Temperature, x_CH4[7], 'lime', linewidth=2, label=labels[7])
    axes3[0, 0].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[0, 0].set_ylabel(r"$x{CH_4}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature);
    lim_max = max(Temperature)
    axes3[0, 0].set_xlim(lim_min, lim_max)
    axes3[0, 0].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(x_CH4[0])
    lim_max = max(x_CH4[7])
    axes3[0, 0].set_ylim(lim_min, lim_max)
    axes3[0, 0].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[0, 0].ticklabel_format(axis='y', style='sci', scilimits=(-5, 1))

    axes3[0, 0].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[0,0].legend(loc='upper left', fontsize=12)

    axes3[1, 0].plot(Temperature, x_CO2[0], 'k', linewidth=2, label=labels[0])
    axes3[1, 0].plot(Temperature, x_CO2[1], 'b', linewidth=2, label=labels[1])
    axes3[1, 0].plot(Temperature, x_CO2[2], 'c', linewidth=2, label=labels[2])
    axes3[1, 0].plot(Temperature, x_CO2[3], 'm', linewidth=2, label=labels[3])
    axes3[1, 0].plot(Temperature, x_CO2[4], 'r', linewidth=2, label=labels[4])
    axes3[1, 0].plot(Temperature, x_CO2[5], 'orange', linewidth=2, label=labels[5])
    axes3[1, 0].plot(Temperature, x_CO2[6], 'brown', linewidth=2, label=labels[6])
    axes3[1, 0].plot(Temperature, x_CO2[7], 'lime', linewidth=2, label=labels[7])
    axes3[1, 0].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[1, 0].set_ylabel(r"$x{CO_2}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes3[1, 0].set_xlim(lim_min, lim_max)
    axes3[1, 0].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(x_CO2[0])
    lim_max = max(x_CO2[7])
    axes3[1, 0].set_ylim(lim_min, lim_max)
    axes3[1, 0].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[1, 0].ticklabel_format(axis='y', style='sci', scilimits=(-4, 4))

    axes3[1, 0].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[1,0].legend(loc='upper left', fontsize=12)

    axes3[2, 0].plot(Temperature, x_H2S[0], 'k', linewidth=2, label=labels[0])
    axes3[2, 0].plot(Temperature, x_H2S[1], 'b', linewidth=2, label=labels[1])
    axes3[2, 0].plot(Temperature, x_H2S[2], 'c', linewidth=2, label=labels[2])
    axes3[2, 0].plot(Temperature, x_H2S[3], 'm', linewidth=2, label=labels[3])
    axes3[2, 0].plot(Temperature, x_H2S[4], 'r', linewidth=2, label=labels[4])
    axes3[2, 0].plot(Temperature, x_H2S[5], 'orange', linewidth=2, label=labels[5])
    axes3[2, 0].plot(Temperature, x_H2S[6], 'brown', linewidth=2, label=labels[6])
    axes3[2, 0].plot(Temperature, x_H2S[7], 'lime', linewidth=2, label=labels[7])
    axes3[2, 0].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[2, 0].set_ylabel(r"$x{H_2S}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes3[2, 0].set_xlim(lim_min, lim_max)
    axes3[2, 0].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(x_H2S[0])
    lim_max = max(x_H2S[7])
    axes3[2, 0].set_ylim(lim_min, lim_max)
    axes3[2, 0].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[2, 0].ticklabel_format(axis='y', style='sci', scilimits=(-4, 4))

    axes3[2, 0].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[2,0].legend(loc='upper left', fontsize=12)

    axes3[3, 0].plot(Temperature, x_H2O[0], 'k', linewidth=2, label=labels[0])
    axes3[3, 0].plot(Temperature, x_H2O[1], 'b', linewidth=2, label=labels[1])
    axes3[3, 0].plot(Temperature, x_H2O[2], 'c', linewidth=2, label=labels[2])
    axes3[3, 0].plot(Temperature, x_H2O[3], 'm', linewidth=2, label=labels[3])
    axes3[3, 0].plot(Temperature, x_H2O[4], 'r', linewidth=2, label=labels[4])
    axes3[3, 0].plot(Temperature, x_H2O[5], 'orange', linewidth=2, label=labels[5])
    axes3[3, 0].plot(Temperature, x_H2O[6], 'brown', linewidth=2, label=labels[6])
    axes3[3, 0].plot(Temperature, x_H2O[7], 'lime', linewidth=2, label=labels[7])
    axes3[3, 0].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[3, 0].set_ylabel(r"$x{H_2O}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes3[3, 0].set_xlim(lim_min, lim_max)
    axes3[3, 0].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(x_H2O[7])
    lim_max = 1
    axes3[3, 0].set_ylim(lim_min, lim_max)
    axes3[3, 0].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[3, 0].ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))

    axes3[3, 0].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[3,0].legend(loc='upper left', fontsize=12)    

    axes3[0, 1].plot(Temperature, y_CH4[0], 'k', linewidth=2, label=labels[0])
    axes3[0, 1].plot(Temperature, y_CH4[1], 'b', linewidth=2, label=labels[1])
    axes3[0, 1].plot(Temperature, y_CH4[2], 'c', linewidth=2, label=labels[2])
    axes3[0, 1].plot(Temperature, y_CH4[3], 'm', linewidth=2, label=labels[3])
    axes3[0, 1].plot(Temperature, y_CH4[4], 'r', linewidth=2, label=labels[4])
    axes3[0, 1].plot(Temperature, y_CH4[5], 'orange', linewidth=2, label=labels[5])
    axes3[0, 1].plot(Temperature, y_CH4[6], 'brown', linewidth=2, label=labels[6])
    axes3[0, 1].plot(Temperature, y_CH4[7], 'lime', linewidth=2, label=labels[7])
    axes3[0, 1].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[0, 1].set_ylabel(r"$y{CH_4}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes3[0, 1].set_xlim(lim_min, lim_max)
    axes3[0, 1].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(y_CH4[0])
    lim_max = max(y_CH4[7])
    axes3[0, 1].set_ylim(lim_min, lim_max)
    axes3[0, 1].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[0, 1].ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))

    axes3[0, 1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[0,1].legend(loc='upper left', fontsize=12)

    axes3[1, 1].plot(Temperature, y_CO2[0], 'k', linewidth=2, label=labels[0])
    axes3[1, 1].plot(Temperature, y_CO2[1], 'b', linewidth=2, label=labels[1])
    axes3[1, 1].plot(Temperature, y_CO2[2], 'c', linewidth=2, label=labels[2])
    axes3[1, 1].plot(Temperature, y_CO2[3], 'm', linewidth=2, label=labels[3])
    axes3[1, 1].plot(Temperature, y_CO2[4], 'r', linewidth=2, label=labels[4])
    axes3[1, 1].plot(Temperature, y_CO2[5], 'orange', linewidth=2, label=labels[5])
    axes3[1, 1].plot(Temperature, y_CO2[6], 'brown', linewidth=2, label=labels[6])
    axes3[1, 1].plot(Temperature, y_CO2[7], 'lime', linewidth=2, label=labels[7])
    axes3[1, 1].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[1, 1].set_ylabel(r"$y{CO_2}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes3[1, 1].set_xlim(lim_min, lim_max)
    axes3[1, 1].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(y_CO2[0])
    lim_max = max(y_CO2[7])
    axes3[1, 1].set_ylim(lim_min, lim_max)
    axes3[1, 1].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[1, 1].ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))

    axes3[1, 1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[1,1].legend(loc='upper left', fontsize=12)

    axes3[2, 1].plot(Temperature, y_H2S[0], 'k', linewidth=2, label=labels[0])
    axes3[2, 1].plot(Temperature, y_H2S[1], 'b', linewidth=2, label=labels[1])
    axes3[2, 1].plot(Temperature, y_H2S[2], 'c', linewidth=2, label=labels[2])
    axes3[2, 1].plot(Temperature, y_H2S[3], 'm', linewidth=2, label=labels[3])
    axes3[2, 1].plot(Temperature, y_H2S[4], 'r', linewidth=2, label=labels[4])
    axes3[2, 1].plot(Temperature, y_H2S[5], 'orange', linewidth=2, label=labels[5])
    axes3[2, 1].plot(Temperature, y_H2S[6], 'brown', linewidth=2, label=labels[6])
    axes3[2, 1].plot(Temperature, y_H2S[7], 'lime', linewidth=2, label=labels[7])
    axes3[2, 1].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[2, 1].set_ylabel(r"$y{H_2S}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature);
    lim_max = max(Temperature)
    axes3[2, 1].set_xlim(lim_min, lim_max)
    axes3[2, 1].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(y_H2S[7])
    lim_max = max(y_H2S[0])
    axes3[2, 1].set_ylim(lim_min, lim_max)
    axes3[2, 1].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[2, 1].ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    axes3[2, 1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[2,1].legend(loc='upper left', fontsize=12)

    axes3[3, 1].plot(Temperature, y_H2O[0], 'k', linewidth=2, label=labels[0])
    axes3[3, 1].plot(Temperature, y_H2O[1], 'b', linewidth=2, label=labels[1])
    axes3[3, 1].plot(Temperature, y_H2O[2], 'c', linewidth=2, label=labels[2])
    axes3[3, 1].plot(Temperature, y_H2O[3], 'm', linewidth=2, label=labels[3])
    axes3[3, 1].plot(Temperature, y_H2O[4], 'r', linewidth=2, label=labels[4])
    axes3[3, 1].plot(Temperature, y_H2O[5], 'orange', linewidth=2, label=labels[5])
    axes3[3, 1].plot(Temperature, y_H2O[6], 'brown', linewidth=2, label=labels[6])
    axes3[3, 1].plot(Temperature, y_H2O[7], 'lime', linewidth=2, label=labels[7])
    axes3[3, 1].set_xlabel(r"$T\,[\mathrm{°C}]$", fontsize=15)
    axes3[3, 1].set_ylabel(r"$y{H_2O}\,[\mathrm{mol/mol}]$", fontsize=15)

    # change by hand every time
    lim_min = min(Temperature)
    lim_max = max(Temperature)
    axes3[3, 1].set_xlim(lim_min, lim_max)
    axes3[3, 1].xaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))

    # change by hand every time
    lim_min = min(y_H2O[7])
    lim_max = max(y_H2O[0])
    axes3[3, 1].set_ylim(lim_min, lim_max)
    axes3[3, 1].yaxis.set_ticks(np.linspace(lim_min, lim_max, 11, endpoint=True))
    axes3[3, 1].ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))

    axes3[3, 1].grid(color='k', alpha=0.8, linestyle='dashed', linewidth=0.8)
    # axes3[3,1].legend(loc='upper left', fontsize=12) 
    axes3[3, 1].legend(loc='upper center', bbox_to_anchor=(-0.14, -0.17), ncol=len(Pressure), fontsize=12)

    path_for_figures = '/Users/fmoretta/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Thermodynamic Optimization ' \
                       'AD/Helmoltz v1/output/VLE figures/'
    filename = path_for_figures + 'x_y_TP_plot.svg'
    fig3.savefig(filename, dpi=800)

    return print('\nsimulation successfull!')
