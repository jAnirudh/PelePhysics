import cantera
from pint import UnitRegistry

ureg = UnitRegistry()
ureg.default_system = "mks"

surf = cantera.Interface("ptcombust.yaml", "Pt_surf")
gas = surf.adjacent["gas"]

surf.TP = 900.0, gas.P

X = dict({key:(1/gas.n_species) for key in gas.species_names})

gas.TPX = 900.0, gas.P, X

surf.advance_coverages(1.0)

gas_npr = gas.net_production_rates
gas_npr *= 1e-3 # convert from kmol/m^3/s to mol/cm^3/s
with open("cantera_data/gas_mole", "w") as logger:
    for idx, spName in enumerate(gas.species_names):
        logger.write(f"wdot_{spName}: {gas_npr[idx]} mol/cm^3/s\n")

gas_npr *= gas.molecular_weights # convert from mol/cm^3 to g/cm^3
with open("cantera_data/gas_mass", "w") as logger:
    for idx, spName in enumerate(gas.species_names):
        logger.write(f"wdot_{spName}: {gas_npr[idx]} g/cm^3/s\n")

surfContri2wdot = surf.get_net_production_rates("gas")
surfWdot = surf.get_net_production_rates("Pt_surf")

surfContri2wdot *= 1e-3 # convert from kmol/m^3/s to mol/cm^3/s
surfWdot *= 0.1 # convert from kmol/m^2/s to mol/cm^2/s

with open("cantera_data/surface_steady", "w") as logger:
    for idx in range(gas.n_species+surf.n_species):
        if not idx < gas.n_species:
            surf_id = idx - gas.n_species
            logger.write(
                f"wdot_{surf.species_names[surf_id]}: {surfWdot[surf_id]}"
                + "mol/cm^2/s\n")
        else:
            logger.write(
                f"wdot_{gas.species_names[idx]}: {surfContri2wdot[idx]}"
                + "mol/cm^3/s\n")

print("Cantera Output is written in cantera_data subdirectory")
