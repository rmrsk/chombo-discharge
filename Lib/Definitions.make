
# These are the names for the various chombo-discharge components. Because we build various components into libraries,
# it makes sense to give the libraries a name through the makefile system. Here, DISCHARGE_LIB is the library name for the
# the library compiled from $DISCHARGE_HOME/Source, GEOMETRY_LIB is the library name for the files compiled from
# $DISCHARGE_HOME/Geometries and so on.
DISCHARGE_LIB  = Discharge
GEOMETRIES_LIB = Geometries

# As a rule we always use EB (embedded boundaries) and MF (multi-fluid) from
# Chombo
USE_EB = TRUE
USE_MF = TRUE

# Chombo libraries needed for building chombo-discharge
LibNames = MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop ParticleTools

# Headers where the chombo-discharge source code is located. This should all folders
# under $DISCHARGE_HOME/Source
SOURCE_DIRS     := $(shell find $(DISCHARGE_HOME)/Source     -type d -print)
GEOMETRIES_DIRS := $(shell find $(DISCHARGE_HOME)/Geometries -type d -print)

# Make variables that hold include flags for various source code.
SOURCE_INCLUDE     := $(foreach dir, $(SOURCE_DIRS),     $(addprefix -I, $(dir)))
GEOMETRIES_INCLUDE := $(foreach dir, $(GEOMETRIES_DIRS), $(addprefix -I, $(dir)))

# Make headers visible
XTRACPPFLAGS += $(SOURCE_INCLUDE) 
XTRACPPFLAGS += $(GEOMETRIES_INCLUDE)

# Make libs visible
XTRALIBFLAGS += $(addprefix -l, $(DISCHARGE_LIB))$(config)
XTRALIBFLAGS += $(addprefix -l, $(GEOMETRIES_LIB))$(config)
XTRALIBFLAGS += -L/$(DISCHARGE_HOME)/Lib
