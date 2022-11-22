z = {
    'CH4': 0.0067,
    'CO2': 0.0022,
    'H2S': 0.0001,
    'H2O': 0.9911
}  # mol/mol

Temperature = [20,  25,  30,  35, 40,  45,  50, 55]  # °C
Pressure = [0.3, 0.4, 0.5, 0.7,  1, 1.2, 1.5, 1.6]  # atm
F = 1000  # kmol/d

dry_z = {
    'CH4': z['CH4']/(sum(z.values())-z['H2O']),
    'CO2': z['CO2']/(sum(z.values())-z['H2O']),
    'H2S': z['H2S']/(sum(z.values())-z['H2O']),
}