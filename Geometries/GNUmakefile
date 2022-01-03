# This file is for building the chombo-discharge source code as a library. The library will be named according to Chombo's naming convection.
include $(DISCHARGE_HOME)/Lib/pre.make

.PHONY: $(GEOMETRIES_LIB) 

$(GEOMETRIES_LIB):

# This is for Chombo -- it will compile the files in src_dirs into a library called GEOMETRIES_LIB
lbase    = $(GEOMETRIES_LIB)
src_dirs = $(GEOMETRIES_DIRS)

# Linker needs to find chombo-discharge source code
include $(DISCHARGE_HOME)/Lib/src.post.make
