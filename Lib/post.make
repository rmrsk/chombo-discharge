# This file should be included at the bottom of every makefile that compiles a chombo-discharge application. If everything works correctly,
# this file will make the chombo-discharge library and header files visible to the preprocessor. Chombo does the rest. 

XTRACPPFLAGS += $(SOURCE_INCLUDE) $(GEOMETRIES_INCLUDE)

XTRALIBFLAGS += $(addprefix -l, $(DISCHARGE_LIB))$(config)
XTRALIBFLAGS += $(addprefix -l, $(GEOMETRY_LIB))$(config)

include $(CHOMBO_HOME)/mk/Make.example
