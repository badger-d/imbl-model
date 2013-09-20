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
  void Process_Ion_Cham_Interacs(std::vector<DetectorHits*> &time_sort_all_hits, vector<DetectorHits*> &time_sort_phot_hits);

  // Method that appends hits to vector.
  void Vectorize_Hits(std::vector<DetectorHits*> &ion_cham_all_hits, std::vector<DetectorHits*> &ion_cham_phot_hits, DetectorHitsCollection* hits_col, std::vector<G4double> &all_hit_times, std::vector<G4double> &phot_hit_times);

  // Method that time sorts the vector of interaction hits.
  void Time_Sort(std::vector<DetectorHits*>	&all_hits, std::vector<G4double> &hit_times, std::vector<DetectorHits*> &time_sort_all_hits);

  // Determine whether or not the photon is a fluorescence or not (primary).
  bool Is_Fluor(const unsigned int track, const unsigned int parent);

  // Determine whether or not the photon is a Rayleigh scatter.
  bool Is_Rayleigh(G4double energy_keV);

  // Determine whether the particle is a photon or electron.
  bool Is_Photon(const G4String particle_ref);

  // Determine whether the electron is from the first primary photon interaction.
  bool Is_Photoelectron(G4String pre_process);

  // Get the energy deposited between the plates.
  bool Energy_Dep_Between_Plates(DetectorHitsCollection* hits_col);

  // Write the primary interaction location to file.
  void Dump_Prim_Photon_Interac_Pos(const G4double x,  const G4double y, const G4double z);

  // Write the energies deposited between the plates to file.
  void Dump_Photon_Interac_Energy(const G4double dump_energy_keV,  const unsigned int tag);

  // Function for collecting the energies deposited in the correct ion chamber layer.
  void Sum_Energies(vector<DetectorHits*> &all_hits, G4double &total_energy_keV, unsigned int layer, G4String volume);

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

  G4int event_id;  // Event ID counter.

  G4double energy_keV_dep_plates_prim; // Total energy deposited  for primary photon interaction between the ionisation plates for the current event.
  G4double energy_keV_dep_plates_scat; // Total energy deposited  for scattered photon interaction between the ionisation plates for the current event.
  G4double energy_keV_dep_plates_fluor; // Total energy deposited  for scattered photon interaction between the ionisation plates for the current event.

  unsigned int sens_layer; // Variable that stores the  ID of the middle layer of the ionisation chamber gas stack .
  G4String sens_vol;       // Variable that stores the string that is linked to the physical volume of the gas stack.

};
#endif
