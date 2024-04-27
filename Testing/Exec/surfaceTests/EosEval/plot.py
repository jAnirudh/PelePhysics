import numpy, pandas
import matplotlib.pyplot as plot

cantera_data = pandas.read_csv("../data/CH4-Pt_EoSTestData.csv")
amrex_data = pandas.read_csv("AMReX_results.csv")

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
