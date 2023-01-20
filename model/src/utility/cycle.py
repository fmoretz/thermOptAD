Temperature = [ 20,  25,  30,  35,  40,  45,  50,  55] # °C
Pressure    = [0.1, 0.2, 0.5, 0.7, 1.0, 1.2, 1.5, 2.0] # atm

i = 0; j = 0

for T in Temperature: 
    
    for P in Pressure:

        print(f'i: {i} -- j: {j}')
        if j == 7:
            j = 0
        else: 
            j = j + 1
    i = i + 1 