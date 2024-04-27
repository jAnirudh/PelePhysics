####################################################################
#            Test data for CEPTR surface evaluations               #
####################################################################
# Author: Anirudh Jonnalagadda                                     #
# email: anirudhj@iisc.ac.in, janirudh@iitb.ac.in                  #
####################################################################

import cantera, numpy, pandas

# basic parameters
p = cantera.one_atm
tinlet = 300.0
tsurface= 900.0
composition = 'CH4:0.095, H2:0.05, O2:0.21, N2:0.78, AR:0.01'

grid = numpy.linspace(0, 1, 129)

def setupImpingingJetCase():
    '''Sets the impinging jet object with the CH4-Pt reaction
    mechanism
    '''
    surf_phase = cantera.Interface("ptcombust.yaml", "Pt_surf")
    surf_phase.TP = tsurface, p
    gas = surf_phase.adjacent["gas"]
    gas.TPX = tinlet, p, composition
    
    surf_phase.advance_coverages(1.0)

    sim = cantera.ImpingingJet(gas=gas, grid=grid,
                               surface=surf_phase)

    sim.set_initial_guess()

    return sim

def writeEosTestData():
    '''Writes csv with test data
    '''
    filename = "CH4-Pt_EoSTestData"
    sim = setupImpingingJetCase()
    # write initial values of grid, velocity, spread-rate,
    # temperature, lambda, and mass-fractions
    sim.write_csv("{}.csv".format(filename))
    # add cp, cv data
    df = pandas.read_csv("{}.csv".format(filename))
    df["cp"] = sim.cp
    df["cv"] = sim.cv
    # overwrite the original csv
    df.to_csv("{}.csv".format(filename))

    with open(filename, "w") as f:
        f.write("\nSurface Cp = {}".format(sim.surface.surface.cp_mole))
        f.write("\nSurface Cv = {}".format(sim.surface.surface.cv_mole))

if __name__ == "__main__":
    writeEosTestData()
