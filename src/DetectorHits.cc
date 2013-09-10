#include "DetectorHits.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<DetectorHits> DetectorHitsAllocator;

//###############################################################################
//Constructor & Destructor

DetectorHits::DetectorHits()
{

}

DetectorHits::~DetectorHits() {
}

//###############################################################################
// Verbose Output

void DetectorHits::Print()
{
    G4cout  << "  Particle: " << particle
	    << "  TrackID: "  << trackID
	    << "  Time: " << time/s
	    << "  ParentID: "  << parentID
	    << "  Position: "  << position/mm
	    << "  Volume: " << volume
	    << "  Delta Position: " << "\t" << delta_position/mm
	    << "  Delta Energy: " << delta_energy/keV
	    << "  Delta momentum: " << "\t" << delta_momentum/keV
	    << "  Process: " << process
    	    << G4endl;
}

//###############################################################################

