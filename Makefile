#-----------------------------------------------------------------------------
# Top level makefile for Bundler
#
# Bundler: Structure from Motion for Unordered Photo Collections
# Version: 0.4 04/03/2010
#    http://phototour.cs.washington.edu/bundler/
#-----------------------------------------------------------------------------
# Copyright (c) 2008-2010 University of Washington and Noah Snavely
# All Rights Reserved.
#-----------------------------------------------------------------------------

ANN_TARGET = linux-g++-shared

OS = $(shell uname -o)
ifeq ($(OS), Cygwin)
ANN_TARGET = win32-g++-shared
endif

default:
# Make libraries
	cd lib/5point; $(MAKE)
	cd lib/ann_1.1_char; $(MAKE) $(ANN_TARGET)
	cd lib/imagelib; $(MAKE)
	cd lib/matrix; $(MAKE)
	cd lib/sba-1.5; $(MAKE)
	cd lib/sfm-driver; $(MAKE)
# Auxiliary libraries
	cd lib/minpack; $(MAKE)
	cd lib/cblas; $(MAKE)
	cd lib/f2c; $(MAKE)
# Main program
	cd src; $(MAKE)


clean:
	cd lib/5point; $(MAKE) clean
	cd lib/ann_1.1_char; $(MAKE) clean
	cd lib/imagelib; $(MAKE) clean
	cd lib/matrix; $(MAKE) clean
	cd lib/sba-1.5; $(MAKE) clean
	cd lib/sfm-driver; $(MAKE) clean
	cd lib/minpack; $(MAKE) clean
	cd lib/cblas; $(MAKE) clean
	cd lib/f2c; $(MAKE) clean
	cd src; $(MAKE) clean
	rm -f bin/bundler bin/KeyMatchFull bin/Bundle2PMVS bin/Bundle2Vis bin/RadialUndistort
	rm -f lib/*.a
