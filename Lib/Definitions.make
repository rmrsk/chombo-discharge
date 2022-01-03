
# These are the names for the various chombo-discharge components. Because we build various components into libraries,
# it makes sense to give the libraries a name through the makefile system. Here, DISCHARGE_LIB is the library name for the
# the library compiled from $DISCHARGE_HOME/Source, GEOMETRY_LIB is the library name for the files compiled from
# $DISCHARGE_HOME/Geometries and so on.
DISCHARGE_LIB = Discharge
GEOMETRY_LIB  = Geometries

# As a rule we always use EB (embedded boundaries) and MF (multi-fluid) from
# Chombo
USE_EB = TRUE
USE_MF = TRUE

# Chombo libraries needed for building chombo-discharge
LibNames = MFElliptic MFTools EBAMRTimeDependent EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools \
	AMRTimeDependent BaseTools BoxTools Workshop ParticleTools


# Headers where the chombo-discharge source code is located. This should all folders
# under $DISCHARGE_HOME/Source
SOURCE_DIRS = $(DISCHARGE_HOME)/Source/AmrMesh \
	$(DISCHARGE_HOME)/Source/ConvectionDiffusionReaction \
	$(DISCHARGE_HOME)/Source/CellTagger \
	$(DISCHARGE_HOME)/Source/Driver \
	$(DISCHARGE_HOME)/Source/Elliptic \
	$(DISCHARGE_HOME)/Source/Electrostatics \
	$(DISCHARGE_HOME)/Source/Geometry \
	$(DISCHARGE_HOME)/Source/ImplicitFunctions \
	$(DISCHARGE_HOME)/Source/Multifluid \
	$(DISCHARGE_HOME)/Source/RadiativeTransfer \
	$(DISCHARGE_HOME)/Source/SigmaSolver \
	$(DISCHARGE_HOME)/Source/Particle \
	$(DISCHARGE_HOME)/Source/Utilities \


# Make a variable which holds all paths to include for source code. 
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

# Headers where the chombo-discharge geometry code is located. This should all folders
# under $DISCHARGE_HOME/Geometries
GEOMETRIES_INCLUDE = -I$(DISCHARGE_HOME)/Geometries/Aerosol \
	-I$(DISCHARGE_HOME)/Geometries/CoaxialCable \
	-I$(DISCHARGE_HOME)/Geometries/DoubleRod \
	-I$(DISCHARGE_HOME)/Geometries/ElectrodeArray \
	-I$(DISCHARGE_HOME)/Geometries/MechanicalShaft \
	-I$(DISCHARGE_HOME)/Geometries/RegularGeometry \
	-I$(DISCHARGE_HOME)/Geometries/RodDielectric \
	-I$(DISCHARGE_HOME)/Geometries/RodPlaneProfile \
	-I$(DISCHARGE_HOME)/Geometries/RoughRod \
	-I$(DISCHARGE_HOME)/Geometries/RoughSphere \
	-I$(DISCHARGE_HOME)/Geometries/Tesselation \
	-I$(DISCHARGE_HOME)/Geometries/Vessel \
	-I$(DISCHARGE_HOME)/Geometries/WireWire \

