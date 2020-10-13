#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (main-level)
#
# !DESCRIPTION: Makefile for the NetCDF Utilities.
#\\
#\\
# !REMARKS:
# To build the program, call "make" with the following syntax:
#
#   make TARGET [ OPTIONAL-FLAGS ]
#
# To display a complete list of options, type "make help".
#
# !REVISION HISTORY: 
#  28 Jul 2015 - R. Yantosca - Initial version, based on GEOS-FP
#EOP
#------------------------------------------------------------------------------
#BOC

#==============================================================================
# Initialization
#==============================================================================

# Define variables
SHELL :=/bin/bash
DIR   :=Code

#==============================================================================
# Makefile targets
#==============================================================================

.PHONY: all lib check clean realclean doc docclean help

all: 
	$(MAKE) -C $(DIR) all

lib:
	$(MAKE) -C $(DIR) lib

check: 
	$(MAKE) -C $(DIR) check

clean:
	$(MAKE) -C $(DIR) clean

realclean:
	$(MAKE) -C $(DIR) realclean

doc:
	$(MAKE) -C $(DIR) doc

docclean:
	$(MAKE) -C $(DIR) docclean

help:
	$(MAKE) -C $(DIR) help

#EOC
