from model import model
from func2 import mole_day
import time

# define operative conditions
Temperature = [20,  25,  30,  35, 40,  45,  50, 55]  # °C
Pressure    = [0.3, 0.5, 1, 1.5]  # atm
species     = ['CH4', 'CO2', 'H2S', 'H2O']  #'H2', 'O2', 'CO', 'NH3', 'ETHANOL', 'METHANOL', 'ETHYLENE']

L = 5 # m first hypothesis
d = 15 # m second hypothesis

# run the model
t1 = time.time()
model(
    Temperature=Temperature,
    Pressure=Pressure,
    species=species
    )

mole_day(
    lks = species[:-1],
    hk  = species[-1],
    Temp_degC = Temperature,
    Pres_atm  = Pressure,
    Height   = L,
    Diameter = d,
    showbool=False
    )
t2 = time.time()

print('\nsimulation successfull!')
print(f'time for execution: {(t2-t1):>2f} seconds')