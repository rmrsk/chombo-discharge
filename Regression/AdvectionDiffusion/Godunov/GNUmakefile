include $(DISCHARGE_HOME)/Lib/pre.make

# Things for the Chombo makefile system. 
ebase    = main
src_dirs = $(DISCHARGE_HOME)/Physics/AdvectionDiffusion

# For building everything. 
all: $(ebase) dependencies

# Build chombo-discharge source code. 
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries

# For building this application -- it needs the chombo-discharge source code. 
$(ebase): dependencies

include $(DISCHARGE_HOME)/Lib/app.post.make
