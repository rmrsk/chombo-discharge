clean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     clean
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries clean

realclean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     realclean
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries realclean

source:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source

geometries: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries

all: source geometries
