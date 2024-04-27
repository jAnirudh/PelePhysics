import subprocess, os

def writeGNUMakeFile(chemistry):
    with open("GNUmakefile", "w") as makefile:
        makefile.write("# AMReX\n")
        makefile.write("DIM       = 3\n")
        makefile.write("PRECISION = DOUBLE\n")
        makefile.write("PROFILE   = FALSE\n")
        makefile.write("VERBOSE   = FALSE\n")
        makefile.write("DEBUG     = FALSE\n")
        makefile.write("\n")
        makefile.write("# Compiler\n")
        makefile.write("COMP     = gnu\n")
        makefile.write("FCOMP    = gfortran\n")
        makefile.write("USE_MPI  = FALSE\n")
        makefile.write("USE_OMP  = FALSE\n")
        makefile.write("USE_CUDA = FALSE\n")
        makefile.write("\n")
        makefile.write("TINY_PROFILE = FALSE\n")
        makefile.write("\n")
        makefile.write("# define the location of the PELE_PHYSICS top directory\n")
        makefile.write("PELE_PHYSICS_HOME    ?= $(shell readlink -f ../../../..)\n")
        makefile.write("AMREX_HOME           ?= $(shell readlink -f ../../../../../amrex)\n")
        makefile.write("\n")
        makefile.write("Eos_Model       = Fuego\n")
        makefile.write("Chemistry_Model = {}\n".format(chemistry))
        makefile.write("Transport_Model = Simple\n")       
        makefile.write("\n")
        makefile.write("Bpack   := ./Make.package\n")
        makefile.write("Blocs   := .\n")
        makefile.write("\n")
        makefile.write("include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics\n")

def run():
    subprocess.call("rm -rf log cantera_data AMReX_data".split(" "))
    os.mkdir("cantera_data")
    os.mkdir("AMReX_data")

    testcases = list([
        "CH4Pt_GASPHASE",
        "chem-CH4-2step-original",
        "chem-CH4-2step"
        ])

    for idx, tcase in enumerate(testcases):

        with open("input", "w") as inputfile:
            inputfile.write(f"testcase = {idx}")

        writeGNUMakeFile(chemistry=tcase)
        subprocess.call("make clean".split(" "))
        subprocess.call("make -j24".split(" "))

        outFile = f"AMReX_data/{tcase}_wdots"

        with open(outFile, "w") as log:
            subprocess.call("./Pele3d.gnu.ex input".split(" "), stdout = log)

        if not tcase == "chem-CH4-2step-original":
            rxn_yaml = f"cantera_test_mechanisms/{tcase}.yaml"
            pyCall = f"python gasPhase_wdot.py -m {rxn_yaml}"
            CanteraOutFile = f"cantera_data/{tcase}_wdots"
            with open(CanteraOutFile, "w") as log:
                subprocess.call(pyCall.split(" "), stdout=log)

        with open("log", "a") as log:
            log.write("*"*(6+len(tcase))+"\n")
            log.write(f"Case: {tcase}\n")
            log.write("*"*(6+len(tcase))+"\n")
            log.write(f"\tAMReX output is written in {outFile}\n")
            log.write(f"\tCantera output is written in {CanteraOutFile}\n")
            log.write("\n")

    subprocess.call("make clean".split(" "))
    subprocess.call("rm GNUmakefile input".split(" "))

    subprocess.call("clear".split(" "))
    subprocess.call("cat log".split(" "))


if __name__ == "__main__":
    run()
