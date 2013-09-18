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
#include "DetectorHits.hh"

class G4Event;
class EventMessenger;
class DetectorConstruction;

class EventAction : public G4UserEventAction
{
public:
  EventAction(PrimaryGeneratorAction*, RunAction*, DetectorConstruction*);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

private:

  // Process the interactions in the ion chamber.
  void Process_Ion_Cham_Interacs(std::vector<DetectorHits*> &all_hits);

  // Method that appends hits to vector.
  void Vectorize_Hits(std::vector<DetectorHits*> &all_hits, DetectorHitsCollection*, std::vector<G4double> &hit_times);

  // Method that time sorts the vector of interaction hits.
  void Time_Sort(std::vector<DetectorHits*>	&all_hits, std::vector<G4double> &hit_times, std::vector<DetectorHits*> &time_sort_all_hits);

  PrimaryGeneratorAction* primary;
  RunAction* run;
  DetectorConstruction* detector;

  DetectorHitsCollection* totHitsCollection;

  G4int ion_cham_layer_hc_id; // Declare hits collection id variable for the ion chamber layers.
  G4int ion_cham_shell_hc_id; // Declare hits collection id variable for the ion chamber shell.
  G4int sample_hc_id;         // Declare hits collection id variable for the sample.
  G4int exp_hall_hc_id;       // Declare hits collection id variable for the experimental hall.

  ofstream* out_file_ptr; // Declare pointer to output file that will be written to.
  G4String out_data;      // Declare string that can contain output file data.

  G4int event_id;  // Event ID conter.

};
#endif
