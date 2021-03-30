# Format of input file
The input file in apothesis is in the form of a .json file. All parameters are assigned with key-value pairs, in the format of {key1:value1, key2:value2}. We use nesting to assign parameters to specific processes or classes. All values should be specified in SI units. Values can be put in any order, but we suggest following the provided template for ease-of-use.

```json
{
    "Lattice":
    {
        "Type": "FCC",
        "dims": [10, 10, 20],
        "species": "Cu"
    }
}
```
In this example, we have defined a copper lattice of type FCC with dimensions [10, 10, 20].

The mandatory keys (specifications) in each input file include
- Lattice
    - Type (FCC, BCC)
    - dims (dimensions of lattice)
    - species
- Temperature
    - In Kelvin
- Pressure
    - In Pa
- Species
    - List of species. See below for example.
- Process
    - Select from "Adsorption", "Desorption", "Diffusion", "Reaction"
    - Each process should further detail important parameters, such as species involved, activation energies
    - See below for descriptions of each process
- config
    - specifies optional configuration parameters such as "debug mode", "writeout rate"
# Lattice
The key of the lattice is "Lattice". 
```json
// @nikos this is repetitive but structured. delete if unnecessary
"Lattice":
{
    "Type": "FCC",
    "dims": [10, 10, 20],
    "species": "Cu"
}
```
Apothesis currently supports lattices of type
- BCC
- FCC 

The dimensions of the lattice specify the size of the lattice in 3 dimensions [x, y, z].

The lattice must must have a species. This is called by assigning a key to the species, usually simply the chemical name. 
# Temperature
Apothesis currently operatures under a constant temperature assumption and does not consider heat effects. The temperature [K] needs to be specified.
```json
Temperature: 400
```
# Pressure
Apothesis operates at a constant pressure assumption. The pressure is specified in Pascals.
```json
Pressure: 101325
```
# Species
All participating reactant species need to be listed here as a nested json structure. For instance, in the case of a CuAMD + H2 -> Cu + H2AMD, we specify:
```json
Species:
{
    "CuAMD": 
    {
        "mw": 130
    },
    "H2":
    {
        "mw": 2
    }
}
```
Each species can optionally provide the molecular weight, in g/mol.
# Process
Each process is specified with a key word. At the moment, the supported processes are:
- Adsorption
- Desorption
- Diffusion
- Reaction

Specific parameters must be defined for each process.
## Adsorption
```json
"Adsorption":
{
    "Species": ["CuAMD", "H2"],
    "Sticking": [0.1, 0.1],
    "MolFraction": [0.67, 0.33],
    "C_tot": 1.03e19,
    "Mass": 32e-3
}
```
The parameters in Adsorption are:
- List of species that can adsorb
- Sticking coefficients [-],
- Energy [kcal/mol]
- mole fractions
- Site density [site/m2]
- Mass [kg/mol]. An arbitrary amount of species can be specified. However, they must always be specified in a list format, even if the case only contains a single species. 

In that case, the adsorption entry would look as follows:
```json
"Adsorption":
{
    "Species": ["CuAMD"],
    "Sticking": [0.1],
    "Energy": [12.0],
    "MolFraction": [0.3],
    "C_tot": 1.0e19,
    "Mass": 32e-3
}
```
The species must be strings, and the values for sticking and mass fraction must be numbers between 0 and 1 (inclusive).
## Desorption
```json
"Desorption":
{
    "Species": ["Cu", "H2"],
    "Energy": [71128, 62220],
    "Frequency":[1e13, 1e13]
}
```
The parameters in desorption are
- Species that can desorb
- Energy of desorption [kcal/mol]
- Diffusion rate when E=E_m=0 [s^-1]
## Diffusion
```json
"Diffusion":
{
    "Species": ["Cu"],
    "E_d": [71128],
    "E_m": [42800],
    "Frequency": [1e13]
}
```
The parameters in diffusion are:
- List of species
- List of desorption energies corresponding to each species [kcal/mole]
- List of energy barrier for each species [kcal/mole]
- Diffusion rate when E=E_m=0 [s^-1]
## Reaction
```json
"Reaction":
{
    "rx1": 
    {
        "Species": ["CuAMD", "H2", "Cu", "HAMD2"],
        "Stoichiometry": [-2, -1, 2, 1],
        "Energy": 66,
        "PreExp": 1e-10
    }
}
```
You are free to define as many reactions as you want. Each reaction should be given a unique name. Here, it is specified as rx1. 

Each listed reaction is defined by:
- List of Species
- Stoichiometric coefficients. Negative values are reactants and positive values are products.
- Energy of reaction [kcal/mol]
- Pre-exponential (Arrhenius) value

It is also posible to only specify the reactants. In that case, the json entry will look as follows:
```json
"Reaction":
{
    "rx1":
    {
        "Species": ["CuAMD", "H2"],
        "Stoichiometry": [-2, -1],
        "Energy": 66,
        "PreExp": 1e-10
    }
}
```
In this case, upon reacting, it is assumed that the products either form the lattice, or diffuse away.
# Common errors
## Comma placement
The json parser is quite strict about comma placement. For instance, the extra comma shown here will result in an error.
```json
"Lattice":
{
    "Type": "FCC",
}
```

## Usage of [] vs individual elements
