# Contains Chombo definitions for makefile system. 
include $(CHOMBO_HOME)/mk/Make.defs

# As a rule we always use EB (embedded boundaries) and MF (multi-fluid) from
# Chombo
USE_EB=TRUE
USE_MF=TRUE

# Chombo libraries needed for building chombo-discharge
LibNames:= MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop ParticleTools

# Headers where the chombo-discharge source code is located. 
SOURCE_INCLUDE = -I./ \
	     -I$(DISCHARGE_HOME)/Source/AmrMesh \
	     -I$(DISCHARGE_HOME)/Source/ConvectionDiffusionReaction \
	     -I$(DISCHARGE_HOME)/Source/CellTagger \
	     -I$(DISCHARGE_HOME)/Source/Driver \
	     -I$(DISCHARGE_HOME)/Source/Elliptic \
	     -I$(DISCHARGE_HOME)/Source/Electrostatics \
	     -I$(DISCHARGE_HOME)/Source/Geometry \
	     -I$(DISCHARGE_HOME)/Source/ImplicitFunctions \
	     -I$(DISCHARGE_HOME)/Source/Multifluid \
	     -I$(DISCHARGE_HOME)/Source/RadiativeTransfer \
	     -I$(DISCHARGE_HOME)/Source/SigmaSolver \
	     -I$(DISCHARGE_HOME)/Source/Particle \
	     -I$(DISCHARGE_HOME)/Source/Utilities \

# # Chombo variable which tells us to look for headers here. 
# XTRACPPFLAGS += $(SOURCE_INCLUDE)
# XTRALIBFLAGS += -lChomboDischarge$(config)
