include $(DISCHARGE_HOME)/Source/Source.make

# Base file containing int main()
ebase    := main
base_dir  = ./
src_dirs  = $(DISCHARGE_HOME)/Geometries/RodDielectric $(DISCHARGE_HOME)/Physics/AdvectionDiffusion

#all: ChomboDischarge $(ebase)
$(ebase): ChomboDischarge

# Build chombo-discharge source code
ChomboDischarge: 
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source ChomboDischarge

# Make sure linker gets all library headers and objects. 
XTRACPPFLAGS += $(SOURCE_INCLUDE)
XTRALIBFLAGS += -lChomboDischarge$(config)

include $(CHOMBO_HOME)/mk/Make.example
