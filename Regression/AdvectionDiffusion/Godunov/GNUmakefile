include $(DISCHARGE_HOME)/Lib/Definitions.make
#include $(DISCHARGE_HOME)/Lib/app.post.make
# Things for the Chombo makefile system. 
ebase    = main
#src_dirs = $(DISCHARGE_HOME)/Physics/AdvectionDiffusion


include $(CHOMBO_HOME)/mk/Make.example

# For building this application -- it needs the chombo-discharge source code. 
$(ebase): dependencies

# Build dependencies if they do not exis. 
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME) lib
	$(MAKE) --directory=$(DISCHARGE_HOME) advectiondiffusion

# Make advection-diffusion headers and library visible. 
XTRACPPFLAGS += $(ADVDIFF_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(ADVDIFF_LIB))$(config)
