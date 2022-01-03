# This file is for building the chombo-discharge source code as a library. The library will be named according to Chombo's naming convection.

include $(DISCHARGE_HOME)/Source/Source.make

.PHONY: ChomboDischarge

lbase := ChomboDischarge

src_dirs=$(DISCHARGE_HOME)/Source/AmrMesh \
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

include $(CHOMBO_HOME)/mk/Make.example
