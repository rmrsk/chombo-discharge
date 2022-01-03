# This file should be included at the bottom of every makefile that needs the chombo-discharge source code. The user needs to ensure
# that the dependency is built before including this file. 

XTRACPPFLAGS += $(SOURCE_INCLUDE)
XTRALIBFLAGS += $(addprefix -l, $(DISCHARGE_LIB))$(config)

include $(CHOMBO_HOME)/mk/Make.example
