# This file is for building the chombo-discharge source code as a library. The library will be named according to Chombo's naming convection.
include $(DISCHARGE_HOME)/Lib/Definitions.make
include $(CHOMBO_HOME)/mk/Make.defs

.PHONY: $(GEOMETRY_LIB)

$(GEOMETRY_LIB): Src

lbase = $(GEOMETRY_LIB)

src_dirs=$(DISCHARGE_HOME)/Geometries/Aerosol\
	$(DISCHARGE_HOME)/Geometries/CoaxialCable	


# Linker needs to find library
XTRACPPFLAGS += $(DISCHARGE_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(DISCHARGE_LIB))$(config)

include $(CHOMBO_HOME)/mk/Make.example
