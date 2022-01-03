include $(DISCHARGE_HOME)/Lib/pre.make

all: $(ebase) dependencies

# Build chombo-discharge dependencies
dependencies: 
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries

# Base file containing int main()
ebase    = main
src_dirs = $(DISCHARGE_HOME)/Physics/AdvectionDiffusion

include $(DISCHARGE_HOME)/Lib/post.make
