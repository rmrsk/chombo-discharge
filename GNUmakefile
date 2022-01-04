clean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries pristine

source:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source lib

geometries: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries lib

advectiondiffusion: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/AdvectionDiffusion lib

brownianwalker: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/BrownianWalker lib

cdrplasma: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/CdrPlasma lib

electrostatics: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Physics/Electrostatics lib

physics: source advectiondiffusion brownianwalker cdrplasma

lib: source geometries 

all: lib advectiondiffusion

