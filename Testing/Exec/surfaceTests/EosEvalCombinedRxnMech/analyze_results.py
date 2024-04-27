import numpy, pandas
import matplotlib.pyplot as plot

cantera_data = pandas.read_csv("../data/CH4-Pt_EoSTestData.csv")
amrex_data = pandas.read_csv("AMReX_GasResults.csv")

x = numpy.linspace(0, 1.0, 129)

plot.plot(x, amrex_data['cp']*1e-4, 'x', label="PelePhysics (X e4)", color="red")
plot.plot(x, cantera_data['cp'], linewidth="3", label="Cantera", color="purple")
plot.legend()
plot.xlabel("x"); plot.ylabel("Cp")
plot.savefig("cp_comparison.png", bbox_inches='tight', dpi=1200)

plot.close()

plot.plot(x, amrex_data['cv']*1e-4, 'x', label="PelePhysics (X e4)", color="red")
plot.plot(x, cantera_data['cv'], linewidth="3", label="Cantera", color="purple")
plot.legend()
plot.xlabel("x"); plot.ylabel("Cv")
plot.savefig("cv_comparison.png", bbox_inches='tight', dpi=1200)

plot.close()


plot.plot(x, amrex_data['density'], 'x', label="PelePhysics", color="red")
plot.plot(x, cantera_data['density'], linewidth="3", label="Cantera", color="purple")
plot.legend()
plot.xlabel("x"); plot.ylabel("Density")
plot.savefig("density_comparison.png", bbox_inches='tight', dpi=1200)

amrex_surface_data = pandas.read_csv("AMReX_SurfaceResults.csv")

with open("../data/CH4-Pt_EoSTestData", "r") as f:
    for line in f.readlines():
        if line.startswith("Surface Cp"):
            cant_surface_cp = line.rstrip()
        if line.startswith("Surface Cv"):
            cant_surface_cv = line.rstrip()

print("\n\nCantera Values (J Kmol^(-1) K^(-1)):\n")
print(cant_surface_cp)
print(cant_surface_cv)

print("\n\nPelePhysics Values (erg mol^(-1) K^(-1)):\n")
print("Cp = {}".format(amrex_surface_data["cp"][128]))
print("Cv = {}".format(amrex_surface_data["cv"][128]))
