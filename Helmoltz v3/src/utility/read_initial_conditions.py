import pandas as pd
from pathlib import Path

cwd = str(Path.cwd())
file_path = cwd + '/'

initial_composition_array = pd.read_excel('./Helmoltz v3/input/input_parameters.xlsx', header=None)

# [column number][row number] - starting from [0][0] = A1
z = {
    'CH4':      float(initial_composition_array[0][3]),
    'CO2':      float(initial_composition_array[1][3]),
    'H2O':      float(initial_composition_array[2][3]),
    'H2S':      float(initial_composition_array[3][3]),
    'H2':       float(initial_composition_array[4][3]),
    'CO':       float(initial_composition_array[5][3]),
    'O2':       float(initial_composition_array[6][3]),
    'NH3':      float(initial_composition_array[7][3]),
    'METHANOL': float(initial_composition_array[8][3]),
    'ETHYLENE': float(initial_composition_array[9][3]),
    'ETHANOL':  float(initial_composition_array[10][3])
}  # mol/mol

F = float(initial_composition_array[2][6])  # flow [kmol/d] - load [kmol]
P = float(initial_composition_array[2][7])  # pressure [atm]
T = float(initial_composition_array[2][8])  # temperature [°C]
L = float(initial_composition_array[2][9])  # height [m]
B = float(initial_composition_array[2][10]) # diameter [m]