import os
import shutil
from decimal import Decimal

stoichiometric_coeff = Decimal(0.9) # yCO = 0.1
stoichiometric_coeff_step = Decimal(0.1) # Increment by 0.1 each time
step = 2

while (stoichiometric_coeff <= Decimal(0.9)):
    file_content = f"""    
    #Build the lattice 
    lattice: SimpleCubic 100 100 10 A 

    #The growing film
    growth: CO2

    #Create steps (alaways in the x-direaction - see above 1st number), number of steps 
    #(e.g. if 10 each step will consist of two atoms), height difference (each step will increase by 20, so 1st step: 10, 2nd: 30, 3rd, 50 etc) 
    #steps: 2 1 

    #Total time in seconds
    time: 1000

    #T in Kelvin
    temperature: 1000

    #P in Pascal
    pressure: 101325

    #Random number initialization
    #random: 123434432

    #Simple s0*f*P/(2*pi*MW*Ctot*kb*T) -> Sticking coefficient [-], f [-], C_tot [sites/m2], MW [kg/mol] 
    #A + * -> A*: simple 0.1 2.0e-4 1.0e+19 0.032 

    #Constant is a constant value in [ML/s]
    CO + * -> CO*: constant {Decimal(stoichiometric_coeff).prec()} #0.389

    #Constant is a constant value in [ML/s]	
    O2 + 2* -> 2O*: constant {Decimal((Decimal(1.0) - Decimal(stoichiometric_coeff)) / Decimal(4.0))} all #0.611
                    
    #A reaction 												
    CO* + O* -> CO2*: constant 1.e+15

    #Example of growth reaction. The 1sr reactant is tranformed to the 1st product, the 2nd reactan to the 2nd product etc. 
    #If #Reactants > #Products then only the first Nth reactants are tranformed to the (N-x) reactants. 
    #Products that have "*" stay on the lattice (adsorbed) 
    #products that do not have "*" are assumed to desorb instantly.  
    #CO* + O* -> CO2* : constant 0.25e+5


    #Time to write in log 
    write: log 0.1

    #Time to write the lattice heights & species
    write: lattice 10

    #Report the coverage of certain species. The time will follow the write in log file
    report: coverage CO* O*
    """

    with open("input.kmc", "w") as f:
        f.write(file_content) # Write the new input.kmc file

    os.system(f"./apothesis") # Run apothesis on the new input.kmc

    # Store the log file into the output file
    shutil.copyfile("Output.log", f"parametric_study/output{step}.txt")

    # Increment the stoichiometric coefficient
    stoichiometric_coeff += stoichiometric_coeff_step
    step += 1