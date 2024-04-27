import subprocess, os

subprocess.call("rm -rf log cantera_data AMReX_data".split(" "))
os.mkdir("cantera_data")
os.mkdir("AMReX_data")

subprocess.call("python net_production_rate.py".split(" "))
subprocess.call("make -j".split(" "))

with open("AMReX_data/solution", "w") as log:
    subprocess.call("./Pele3d.gnu.ex".split(" "), stdout=log)

subprocess.call("make clean".split(" "))
