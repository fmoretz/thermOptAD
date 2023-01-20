z = {
    'CH4': 0.0067,
    'CO2': 0.0022,
    'H2S': 0.0001,
    'H2O': 0.9911,
    'H2': 0,
    'O2': 0,
    'CO': 0,
    'NH3': 0,
    'ETHANOL': 0,
    'METHANOL': 0,
    'ETHYLENE': 0
}  # mol/mol

spec_sel = ['CH4', 'ETHYLENE', 'O2', 'METHANOL']

reduced_z = {key: z[key] for key in spec_sel}

print(reduced_z)