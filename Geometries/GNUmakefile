# Contains various definitions of which Chombo libraries we need, and where the various chombo-discharge source files are located
include $(DISCHARGE_HOME)/Lib/Definitions.make

# This is for Chombo -- it will compile the files in src_dirs into a library called $DISCHARGE_LIB. Chombo's makefiles will expand
# all files found in src_dirs into $objects. 
src_dirs = $(GEOMETRIES_DIRS)
include $(CHOMBO_HOME)/mk/Make.example

# Name of the library
SrcLib = $(DISCHARGE_HOME)/Lib/lib$(GEOMETRIES_LIB)$(config).a

# For making $DISCHARGE_HOME/Geometries into a library in $DISCHARGE_HOME/Lib
lib: $(objects)
	if test -f $(SrcLib); then $(RM) $(SrcLib); fi
	ar -cr $(SrcLib) o/$(config)/*.o

all: lib

pristine:
	$(RM) -r [odpf]/$(config)
	if test -f $(SrcLib); then $(RM) $(SrcLib); fi 

.PHONY: all lib pristine
