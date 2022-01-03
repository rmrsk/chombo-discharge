# This file is for building the chombo-discharge source code as a library. The library will be named according to Chombo's naming convection.
include $(DISCHARGE_HOME)/Lib/pre.make

.PHONY: $(DISCHARGE_LIB) 

$(DISCHARGE_LIB):

# This is for Chombo -- it will compile the files in src_dirs into a library called $DISCHARGE_LIB
lbase    = $(DISCHARGE_LIB)
src_dirs = $(SOURCE_DIRS)

# Chombo does it's thing. 
include $(CHOMBO_HOME)/mk/Make.example
