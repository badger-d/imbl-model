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

  // Determine whether or not the photon is a primary or not (fluorescence).
  bool Is_Primary(const unsigned int track, const unsigned int parent);

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
  void Dump_Photon_Interac_Energy(const G4double dump_energy_keV,  const unsigned int output_flag);

  // Function for collecting the energies that are deposited in the correct ion chamber layer.
  G4double Sum_Energies_In_Hit_Vol(std::vector<DetectorHits*> &all_hits, unsigned int sensitive_layer, G4String sensitive_volume);

  // Function for collecting all other energies are not deposited in the sensitive ion chamber layer.
  G4double Sum_Energies_Not_In_Hit_Vol(std::vector<DetectorHits*> &all_hits, unsigned int sensitive_layer, G4String sensitive_volume, G4String shell_volume);

  // Check if all of the electron deposits were in the gas.
  bool All_Electron_Interacs_In_Sens(std::vector<DetectorHits*> &all_hits);

  // Get the progeny of a photon interaction for photoelectric absorption.
  void Get_Photon_Progeny_Phot(vector<DetectorHits*> &all_hits, unsigned int init_interac_index, vector<DetectorHits*> &daughter_tracks);

  // Get the progeny of a photon interaction for Compton scatter.
  void Get_Photon_Progeny_Scat(vector<DetectorHits*> &all_hits, unsigned int init_interac_index, vector<DetectorHits*> &daughter_tracks);

  // Get the progeny of a fluorescence photon interactions.
  void Get_Photon_Progeny_Fluor(vector<DetectorHits*> &all_hits, vector<unsigned int> &indices, vector<DetectorHits*> &daughter_tracks);

  // Add the track ID to the vector.
  void Add_Track_ID(unsigned int track, vector<unsigned int> &parent_vector);

  // Determine whether the current interaction has a parent and is therefore progeny.
  bool Has_Relevant_Parent(unsigned int parent, vector<unsigned int> parent_vector);

  // Get the first primary photon interaction.
  unsigned int Get_First_Photon_Index(vector<DetectorHits*> &all_hits);

  // Determine whether or not a scatter has occurred.
  bool Has_Scat(std::vector<DetectorHits*> &phot_hits);

  // Determine whether or not a fluorescence photon is present.
  bool Has_Fluor(vector<DetectorHits*> &phot_hits);

  // Get vector of interaction indices for fluorescence interactions.
  void Get_Fluor_Photon_Indices(vector<DetectorHits*> &all_hits, vector<unsigned int> &indices);

  // Determine whether part of the scatter event has left the sensitive volume.
  bool Has_Left_Sens(vector<DetectorHits*> &phot_hits, unsigned int layer, G4String volume);

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

  G4double energy_keV_dep_plates_prim_vola;  // Total energy deposited for primary photon interaction between the ionisation plates for the current event.
  G4double energy_keV_dep_plates_prim_volb;  // Total energy deposited for primary photon interaction in regions other than between the ionisation plates for the current event.
  G4double energy_keV_dep_plates_scat_vola;  // Total energy deposited for scattered photon interaction between the ionisation plates for the current event.
  G4double energy_keV_dep_plates_fluor_vola; // Total energy deposited for scattered photon interaction between the ionisation plates for the current event.
  G4double energy_keV_dep_plates_eloss_vola; // Total energy deposited for photons whose progeny reached the collecting electrode - volume A.
  G4double energy_keV_dep_plates_eloss_volb; // Total energy deposited for photons whose progeny reached the collecting electrode - volume B.

  unsigned int sens_layer; // Variable that stores the  ID of the middle layer of the ionisation chamber gas stack .
  G4String sens_vol;       // Variable that stores the string that is linked to the physical volume of the gas stack.
  G4String shell_vol;      // Variable that stores the string that is linked to the physical volume of the ion chamber shell.
};
#endif
