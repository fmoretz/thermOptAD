import matplotlib.pyplot as plt
from sympy import sympify
from func2 import *
import numpy as np
from initial import *
from pathlib import Path

def model(Temperature, Pressure, species):
    plt.style.use(['science','no-latex'])
    alpha0 = z['CH4'] + z['CO2'] + z['H2S']

    alpha_TP = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    x_CH4    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    y_CH4    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    x_CO2    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    y_CO2    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    x_H2S    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    y_H2S    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    x_H2O    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    y_H2O    = np.empty(shape=(len(Pressure), len(Temperature)))  # --
    Vap_TP   = np.empty(shape=(len(Pressure), len(Temperature)))  # kmol/d
    Liq_TP   = np.empty(shape=(len(Pressure), len(Temperature)))  # kmol/d

    i = 0  # index for cycling in T for-loop
    j = 0  # index for cycling in P for-loop

    for T in Temperature:
        T = T + 273  # K

        for P in Pressure:
            P = P * 101325  # Pa

            # initialization
            Ar_V       = {key: 0 for key in species}  # --
            entropy_V  = {key: 0 for key in species}  # J/mol/K
            enthalpy_V = {key: 0 for key in species}  # J/mol/K
            cp_V       = {key: 0 for key in species}  # J/mol/K
            cv_V       = {key: 0 for key in species}  # J/mol/K
            cubicExp_V = {key: 0 for key in species}  # --
            fuga_V     = {key: 0 for key in species}  # --
            Ar_L       = {key: 0 for key in species}  # --
            entropy_L  = {key: 0 for key in species}  # J/mol/K
            enthalpy_L = {key: 0 for key in species}  # J/mol/K
            cp_L       = {key: 0 for key in species}  # J/mol/K
            cv_L       = {key: 0 for key in species}  # J/mol/K
            cubicExp_L = {key: 0 for key in species}  # --
            fuga_L     = {key: 0 for key in species}  # --
            Psat       = {key: 0 for key in species}  # Pa
            dens_L     = {key: 0 for key in species}  # kg * m^(-3)
            dens_V     = {key: 0 for key in species}  # kg * m^(-3)
            x          = {key: 0 for key in species}  # --
            y          = {key: 0 for key in species}  # --
            Lx         = {key: 0 for key in species}  # kmol/d
            Vy         = {key: 0 for key in species}  # kmol/d

            for specie in species:

                if specie == 'CO':
                    Psat[specie] = Antoine(Ant[specie], T-273.15) * 133.322  # mmHg to Pa

                if specie == 'CH4':
                    Psat[specie] = Henry(kH[specie], T) * 1e3 # kPa to Pa

                if specie == 'CO2':
                    Psat[specie] = Henry(kH[specie], T) * 1e3 # kPa to Pa

                else:
                    Psat[specie] = Antoine(Ant[specie], T) * 101325 # kPa to Pa

                dens_L[specie]     = density(d_liq[specie], T, Psat[specie], specie)
                dens_V[specie]     = density(d_vap[specie], T, P, specie)
                Ar_L[specie]       = sympify(VLProperties(0, T, dens_L[specie], specie))
                entropy_L[specie]  = sympify(VLProperties(1, T, dens_L[specie], specie))
                enthalpy_L[specie] = sympify(VLProperties(2, T, dens_L[specie], specie))
                cp_L[specie]       = sympify(VLProperties(3, T, dens_L[specie], specie))
                cv_L[specie]       = sympify(VLProperties(4, T, dens_L[specie], specie))
                cubicExp_L[specie] = sympify(VLProperties(5, T, dens_L[specie], specie))
                fuga_L[specie]     = sympify(VLProperties(6, T, dens_L[specie], specie))
                Ar_V[specie]       = sympify(VLProperties(0, T, dens_V[specie], specie))
                entropy_V[specie]  = sympify(VLProperties(1, T, dens_V[specie], specie))
                enthalpy_V[specie] = sympify(VLProperties(2, T, dens_V[specie], specie))
                cp_V[specie]       = sympify(VLProperties(3, T, dens_V[specie], specie))
                cv_V[specie]       = sympify(VLProperties(4, T, dens_V[specie], specie))
                cubicExp_V[specie] = sympify(VLProperties(5, T, dens_V[specie], specie))
                fuga_V[specie]     = sympify(VLProperties(6, T, dens_V[specie], specie))

            epsilon = 1
            x_guess = {key: z[key] for key in species}
            while epsilon > 1e-3:
                fuga_m = {key: 0 for key in species}  # --
                for specie in species:
                    fuga_m[specie] = fug_mix_helmoltz(specie, T, P, x_guess)
                k = {
                    key: Psat[key] * fuga_L[key] / (P * fuga_m[key]) for key in species
                } # Pa/Pa

                # Rachford-Rice Solution
                reduced_z = {key: z[key] for key in species}
                reduced_k = {key: k[key] for key in species}

                alpha = fsolve(
                    func = RR,
                    x0   = alpha0,
                    args = (
                        list(reduced_z.values()),
                        list(reduced_k.values())
                    )
                )

                alpha = alpha[0]
                V = alpha * F
                L = F - V

                for specie in species:
                    x[specie]  = z[specie] / (1 + alpha * (k[specie] - 1))
                    y[specie]  = x[specie] * k[specie]
                    Vy[specie] = V * y[specie]
                    Lx[specie] = L * x[specie]

                eps = {key: np.abs(x_guess[key] - y[key]) for key in species}
                epsilon = np.mean(list(eps.values()))
                print(f"T:{T-273}°C | P:{P/101325} atm | eps: {epsilon}")
                print(fuga_m)
                x_guess = y

            # array allocation of pressure(j)-temperature(i) dependent variables
            alpha_TP[j][i] = alpha
            x_CH4[j][i]  = x['CH4']
            y_CH4[j][i]  = y['CH4']
            x_CO2[j][i]  = x['CO2']
            y_CO2[j][i]  = y['CO2']
            x_H2S[j][i]  = x['H2S']
            y_H2S[j][i]  = y['H2S']
            x_H2O[j][i]  = x['H2O']
            y_H2O[j][i]  = y['H2O']
            Vap_TP[j][i] = V
            Liq_TP[j][i] = L


            # Saving thermodynamic properties in a txt file
            P_str = str(P / 101325)
            P_str = P_str.replace('.', '')

            filename = f'results_{int(T - 273)}C_' + P_str + 'atm'
            # cd = os.chdir('')
            path = '/output/thermo parameters/'

            filename = path + filename + '.txt'

            # retrieve file path
            # `cwd`: current directory
            cwd = Path.cwd()
            file = str(cwd) + filename

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

            # Saving compositions in a txt file
            P_str = str(P / 101325)
            P_str = P_str.replace('.', '')

            filename = f'xy_{int(T - 273)}C_' + P_str + 'atm'
            # cd = os.chdir('')
            path = '/output/compositions/txt/'

            filename = path + filename + '.txt'

            file = str(cwd) + filename

            with open(file, 'w') as f:
                f.write("xCH4\txCO2\txH2S\txH2O\tyCH4\tyCO2\tyH2S\tyH2O\n")
                f.write(f"{x_CH4[j][i]}\t{x_CO2[j][i]}\t{x_H2S[j][i]}\t{x_H2O[j][i]}\
                    \t{y_CH4[j][i]}\t{y_CO2[j][i]}\t{y_H2S[j][i]}\t{y_H2O[j][i]}\n")

            if j == len(Pressure)-1:
                j = 0
            else:
                j = j + 1
        i = i + 1

    labels = [
        'P = ' + str(Pressure[0]) + ' atm',
        'P = ' + str(Pressure[1]) + ' atm',
        'P = ' + str(Pressure[2]) + ' atm',
        'P = ' + str(Pressure[3]) + ' atm'
    ]

    # visualization settings
    tick_num = 10
    dpi = 250
    x_lim_min = min(Temperature)
    x_lim_max = max(Temperature)
    Temp_label = r"$\mathrm{T\,[°C]}$"
    lw = 0.75
    ms = 2

    # alpha plot
    fig, axes = plt.subplots(1, 1, dpi=dpi)
    axes.plot(Temperature, alpha_TP[0], marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, alpha_TP[1], marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, alpha_TP[2], marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, alpha_TP[3], marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{\gamma\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(0.003, 0.02)
    plt.ticklabel_format(axis='y', style='sci', scilimits=[-2, 0])
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    #axes.legend(loc='best')

    path_for_figures = '/output/VLE figures/'
    file = 'alpha_TP_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    # V-L plot
    fig, axes = plt.subplots(1, 1, dpi=dpi)
    axes.plot(Temperature, Vap_TP[0], marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, Vap_TP[1], marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, Vap_TP[2], marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, Vap_TP[3], marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{Vapour\,flow\,[kmol/d]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(Vap_TP[len(Pressure)-1])/1.02, max(Vap_TP[0]))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'V_L_TP_plot_1.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi=dpi)
    axes.plot(Temperature, Liq_TP[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, Liq_TP[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, Liq_TP[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, Liq_TP[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{Liquid\,flow\,[kmol/d]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(Liq_TP[0]) - 2, F)
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    #axes.legend(loc='best')

    file = 'V_L_TP_plot_2.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    # x-y plot
    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, x_CH4[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, x_CH4[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, x_CH4[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, x_CH4[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{x_{CH_4}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(x_CH4[0]), max(x_CH4[len(Pressure)-1]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'x_y_TP_xCH4_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, x_CO2[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, x_CO2[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, x_CO2[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, x_CO2[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{x_{CO_2}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(x_CO2[0]), max(x_CO2[len(Pressure)-1]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-2, 0))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    #axes.legend(loc='best')

    file = 'x_y_TP_xCO2_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, x_H2S[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, x_H2S[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, x_H2S[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, x_H2S[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{x_{H_2S}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(x_H2S[0]), max(x_H2S[len(Pressure)-1]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-4, 4))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'x_y_TP_xH2S_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, x_H2O[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, x_H2O[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, x_H2O[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, x_H2O[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{x_{H_2O}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(x_H2O[len(Pressure)-1]), 1)
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'x_y_TP_xH2O_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, y_CH4[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, y_CH4[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, y_CH4[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, y_CH4[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{y_{CH_4}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(y_CH4[0]), max(y_CH4[len(Pressure)-1]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'x_y_TP_yCH4_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, y_CO2[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, y_CO2[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, y_CO2[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, y_CO2[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{y_{CO_2}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(y_CO2[0]), max(y_CO2[len(Pressure)-1]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'x_y_TP_yCO2_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, y_H2S[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, y_H2S[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, y_H2S[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, y_H2S[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{y_{H_2S}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(y_H2S[len(Pressure)-1]), max(y_H2S[0]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)

    file = 'x_y_TP_yH2S_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)

    fig, axes = plt.subplots(1, 1, dpi = dpi)
    axes.plot(Temperature, y_H2O[0],  marker='o', markersize=ms, linewidth=lw, label=labels[0])
    axes.plot(Temperature, y_H2O[1],  marker='^', markersize=ms, linewidth=lw, label=labels[1])
    axes.plot(Temperature, y_H2O[2],  marker='s', markersize=ms, linewidth=lw, label=labels[2])
    axes.plot(Temperature, y_H2O[3],  marker='D', markersize=ms, linewidth=lw, label=labels[3])
    axes.set_xlabel(Temp_label)
    axes.set_ylabel(r"$\mathrm{y_{H_2O}\,[--]}$")
    axes.set_xlim(x_lim_min, x_lim_max)
    axes.set_ylim(min(y_H2O[len(Pressure)-1]), max(y_H2O[0]))
    axes.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    axes.grid(color='k', alpha=0.2, linestyle='dashed', linewidth=0.5)
    #axes.legend(loc='best')

    file = 'x_y_TP_yH2O_plot.svg'
    filename = str(cwd) + path_for_figures + file
    fig.savefig(filename, dpi=dpi)