from model import model
from func2 import mole_day
import time
from initial import T, P

# define operative conditions
Temperature = T # °C
Pressure    = P # atm
species     = ['CH4', 'CO2', 'H2S', 'H2O']  #'H2', 'O2', 'CO', 'NH3', 'ETHANOL', 'METHANOL', 'ETHYLENE']

# run the model
t1 = time.time()
model(
    Temperature=Temperature,
    Pressure=Pressure,
    species=species
    )

t2 = time.time()

print('\nsimulation successfull!')
print(f'time for execution: {(t2-t1):>2f} seconds')