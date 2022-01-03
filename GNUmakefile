clean:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source     pristine
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries pristine

source:
	$(MAKE) --directory=$(DISCHARGE_HOME)/Source lib

geometries: source
	$(MAKE) --directory=$(DISCHARGE_HOME)/Geometries lib

lib: source geometries

