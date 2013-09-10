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

private:

  PrimaryGeneratorAction* primary;
  RunAction* run;
  DetectorConstruction* detector;

  G4int ion_cham_hc_id; // Declare hits collection id variable for the ion chamber.
  G4int sample_hc_id;   // Declare hits collection id variable for the sample.
  G4int exp_hall_hc_id; // Declare hits collection id variable for the experimental hall.

  ofstream* out_file_ptr; // Declare pointer to output file that will be written to.
  G4String out_data;      // Declare string that can contain output file data.

  G4int event_id;  // Event ID conter.

};
#endif
