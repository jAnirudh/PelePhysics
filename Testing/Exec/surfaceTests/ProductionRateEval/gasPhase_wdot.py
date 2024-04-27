from pint import UnitRegistry
import numpy, cantera, argparse, os

cantera.use_legacy_rate_constants(False)
ureg = UnitRegistry()
ureg.default_system = 'cgs'

# Cantera Units for molar density
rho_cant_units = ureg.kilomole/ureg.m**3
# Base (cgs) Units for molar density
rho_base_units = (1*rho_cant_units).to_base_units().units
# Conversion Factor for molar density
rho_ConvFac = 1/(1*rho_cant_units).to_base_units().m

# Cantera units for progress rate
q_cant_units = ureg.kilomole/ureg.m**3/ureg.s
# Base (CGS) units for progress rate
q_base_units = (1*q_cant_units).to_base_units().units
# Conversion Factor for progress rate
q_ConvFac = (1*q_cant_units).to_base_units().m

# Cantera units for molecular weights
mw_cant_units = ureg.kilogram/ureg.kilomole
# Base (CGS) units for molecular weights
mw_base_units = (1*mw_cant_units).to_base_units().units
# Conversion Factor for molecular weights
mw_ConvFac = 1/(1*mw_cant_units).to_base_units().m

# empty lists for storing forward/reverse rate constant
kf_cant_units = list() ; kr_cant_units = list()
kf_base_units = list() ; kr_base_units = list()
# empty lists for storing creation/destruction rates of progress
qf_cant_units = list() ; qr_cant_units = list()
qf_base_units = list() ; qr_base_units = list()

def get_unit_info(reaction):
    # Cantera Units for pre-exponential constant
    A_cant_units_str = str(reaction.rate_coeff_units).split("at")[0][7:-2]
    A_cant_units     = ureg.parse_expression(A_cant_units_str).units
    #CGS units (base)
    A_base_units   = (1*  A_cant_units).to_base_units().units
    # conversion factors
    A_ConvFac   = 1/(1*  A_cant_units).to_base_units().m
    return A_cant_units, A_base_units, A_ConvFac

def run(mechanism):
    gas = cantera.Solution(mechanism)
    if mechanism.endswith("GASPHASE.yaml"):
        X = dict({key:(1/gas.n_species) for key in gas.species_names})
        T = gas.T # 300K
    else:
        X = "CH4: 0.5, O2: 0.5"
        T = 900 # 900K

    gas.TPX = T, gas.P, X # p =1atm

    # molar net production rates
    molar_npr_cant_units = list([pr*q_cant_units for pr in gas.net_production_rates])
    molar_npr_base_units = list([pr*q_ConvFac*q_base_units for pr in gas.net_production_rates])

    # net production rates
    npr_cant_units = list([m_npr*mw*mw_cant_units for m_npr, mw in zip(molar_npr_cant_units, gas.molecular_weights)])
    npr_base_units = list([m_npr*mw*mw_base_units for m_npr, mw in zip(molar_npr_base_units, gas.molecular_weights)])

    for idx, sp_name in enumerate(gas.species_names):
        print(f"wdot_{sp_name}: {npr_base_units[idx].m} {npr_base_units[idx].units}")

if __name__ == "__main__":
    if not os.path.exists(os.getcwd()+"/cantera_test_mechanisms/"):
        raise(FileNotFoundError,
              '"cantera_test_mechanisms" folder not found in CWD')

    parser = argparse.ArgumentParser(description='Compute Cantera Production rates')
    group = parser.add_mutually_exclusive_group()
    parser.add_argument('-m', '--mechanismFile', type=str,
                        metavar='', required=False,
                        help='Reaction Mechanism',
                        default='cantera_test_mechanisms/CH4-2step-mod.yaml')
    args = parser.parse_args()

    if not os.path.exists(args.mechanismFile):
        raise FileNotFoundError('"{}" not found.'.format(args.mechanismFile))

    run(mechanism = args.mechanismFile)
