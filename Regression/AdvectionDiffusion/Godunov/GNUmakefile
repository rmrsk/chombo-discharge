include $(CHOMBO_HOME)/mk/Make.defs

# Base file containing int main()
ebase    = main
src_dirs = $(DISCHARGE_HOME)/Geometries/RodDielectric $(DISCHARGE_HOME)/Physics/AdvectionDiffusion

# # Build this application. 
$(ebase): source

# # Build chombo-discharge source code
source: 
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source

# Make sure linker gets all library headers and objects. 
XTRACPPFLAGS += $(SOURCE_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(DISCHARGE_LIB))$(config)

include $(CHOMBO_HOME)/mk/Make.example
