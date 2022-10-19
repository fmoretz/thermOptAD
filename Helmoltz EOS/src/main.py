from model import model

# define operative conditions
Temperature = [ 20,  25,  30,  35, 40,  45,  50, 55] # °C
Pressure    = [0.3, 0.4, 0.5, 0.7,  1, 1.2, 1.5,  2] # atm

# run the model
model(
    Temperature = Temperature,
    Pressure    = Pressure
    )
