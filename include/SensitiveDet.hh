#ifndef SensitiveDet_h
#define SensitiveDet_h 1

#include "G4VSensitiveDetector.hh"
#include "DetectorHits.hh"
#include "globals.hh"
#include <fstream>
#include <vector>

using namespace std;
class G4Step;
class G4HCofThisEvent;
class DetectorConstruction;

//----------------------------------------------------------------------------

class SensitiveDet : public G4VSensitiveDetector
{
  public:
      SensitiveDet(G4String);
     virtual ~SensitiveDet();

     virtual void Initialize(G4HCofThisEvent*HCE);
     virtual G4bool ProcessHits(G4Step*aStep, G4TouchableHistory*ROhist);
     virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:

      void SetVolume(G4String val){volume = val;};
      G4String GetVolume(){return volume;};

      G4String volume;
      DetectorHitsCollection* hitsCollection;
      ofstream* foutPositions;  //
      G4int particles;                      //
      DetectorConstruction* Detector;
      G4int HCID;

      G4int parentID;
      G4String particle;

};

#endif
