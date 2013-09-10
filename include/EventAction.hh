//filename : EventAction.hh
//----------------------------------------------------------------------------
#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <fstream>
#include "DetectorHits.hh"
#include "PrimaryGeneratorAction.hh"
#include <vector>
#include <algorithm>
#include "RunAction.hh"
#include "DetectorConstruction.hh"

class G4Event;
class EventMessenger;
class DetectorConstruction;

using namespace std;

class EventAction : public G4UserEventAction
{
public:
  EventAction(PrimaryGeneratorAction*, RunAction*, DetectorConstruction*);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  //messenger functions

private:

  PrimaryGeneratorAction* primary;
  RunAction* run;
  DetectorConstruction* detector;

  G4int detectorHC_ID;
  G4int sampleHC_ID;
  G4int hallHC_ID;

  ofstream* outFile_ptr;

  G4String outData;

};
#endif
