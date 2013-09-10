# --------------------------------------------------------------
# GNUmakefile for bones.
# --------------------------------------------------------------

CPPVERBOSE := true

name := silicon
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib bin

#include $(G4INSTALL)/config/binmake.gmk
include $(G4INSTALL)/share/Geant4-9.5.0/geant4make/geant4make.sh
CPPFLAGS += -Wno-deprecated -g3 -O3

# When running on cluster, uncomment these lines.
#CPPFLAGS += -Wno-deprecated -g3 -O3 -I/home/matthew/boost/include
#LOADLIBS += -L/home/matthew/boost/lib

