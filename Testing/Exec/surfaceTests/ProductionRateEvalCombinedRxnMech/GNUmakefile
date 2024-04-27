# AMReX
DIM       = 3
PRECISION = DOUBLE
PROFILE   = FALSE
VERBOSE   = FALSE
DEBUG     = FALSE

# Compiler
COMP     = gnu
FCOMP    = gfortran
USE_MPI  = FALSE
USE_OMP  = FALSE
USE_CUDA = FALSE
USE_HIP  = FALSE

HETEROGENEOUS_CHEMISTRY = TRUE
STEADY_STATE = TRUE
TINY_PROFILE = FALSE

# define the location of the PELE_PHYSICS top directory
PELE_PHYSICS_HOME    ?= $(shell readlink -f ../../../..)
AMREX_HOME           ?= $(shell readlink -f ../../../../../amrex)

Eos_Model       = Fuego
Chemistry_Model = CH4Pt
Transport_Model = Constant
Surface_Model = Ideal


Bpack   := ./Make.package
Blocs   := .

include $(PELE_PHYSICS_HOME)/Testing/Exec/Make.PelePhysics
