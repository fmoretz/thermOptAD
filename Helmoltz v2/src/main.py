from model import model

# define operative conditions
Temperature = [20,  25,  30,  35, 40,  45,  50, 55]  # °C
Pressure = [0.3, 0.4, 0.5, 0.7,  1, 1.2, 1.5, 1.6]  # atm
species = ['CH4', 'CO2', 'H2S', 'H2O'] #'H2', 'O2', 'CO', 'NH3', 'ETHANOL', 'METHANOL', 'ETHYLENE']
# run the model
model(
    Temperature=Temperature,
    Pressure=Pressure,
    species=species
    )
