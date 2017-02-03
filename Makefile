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
PREFIX = $(DESTDIR)/usr
BINDIR = $(PREFIX)/bin
LIBDIR = $(PREFIX)/lib

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

install:
	cd bin;
	mkdir -p $(BINDIR)
	mkdir -p $(LIBDIR)

# Install shared library
	install -Dm644 bin/libANN_char.so $(LIBDIR)/libANN_char.so

# Install binaries
	install -Dm755 bin/Bundle2Ply $(BINDIR)/Bundle2Ply
	install -Dm755 bin/Bundle2PMVS $(BINDIR)/Bundle2PMVS
	install -Dm755 bin/Bundle2Vis $(BINDIR)/Bundle2Vis
	install -Dm755 bin/bundler $(BINDIR)/bundler_sfm
	install -Dm755 bin/extract_focal.pl $(BINDIR)/extract_focal.pl
	install -Dm755 bin/FisheyeUndistort $(BINDIR)/FisheyeUndistort
	install -Dm755 bin/KeyMatchFull $(BINDIR)/KeyMatchFull
	install -Dm755 bin/RadialUndistort $(BINDIR)/RadialUndistort
	install -Dm755 bin/ToSift.sh $(BINDIR)/ToSift.sh
	install -Dm755 bin/ToSiftList.sh $(BINDIR)/ToSiftList.sh

uninstall:
	# Uninstall shared library
	rm $(LIBDIR)/libANN_char.so

	# Uninstall binaries
	rm $(BINDIR)/Bundle2Ply
	rm $(BINDIR)/Bundle2PMVS
	rm $(BINDIR)/Bundle2Vis
	rm $(BINDIR)/bundler_sfm
	rm $(BINDIR)/extract_focal.pl
	rm $(BINDIR)/FisheyeUndistort
	rm $(BINDIR)/KeyMatchFull
	rm $(BINDIR)/RadialUndistort
	rm $(BINDIR)/ToSift.sh
	rm $(BINDIR)/ToSiftList.sh
	
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
