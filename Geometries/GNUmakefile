# This file is for building the chombo-discharge source code as a library. The library will be named according to Chombo's naming convection.
include $(DISCHARGE_HOME)/Lib/Definitions.make
include $(CHOMBO_HOME)/mk/Make.defs

lbase = $(GEOMETRIES_LIB)

src_dirs=$(GEOMETRIES_DIRS)

# Linker needs to find library
XTRACPPFLAGS += $(SOURCE_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(DISCHARGE_LIB))$(config)

include $(CHOMBO_HOME)/mk/Make.example
