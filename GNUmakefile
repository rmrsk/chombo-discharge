clean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     realclean
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries realclean

source:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source $(DISCHARGE_LIB)

geometries: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries $(GEOMETRY_LIB)

all: source geometries
