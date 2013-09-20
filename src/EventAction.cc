#include "G4Event.hh"
#include "EventAction.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

using namespace std;

EventAction::EventAction(PrimaryGeneratorAction* PGA, RunAction* RA, DetectorConstruction* DC)
{
   primary = PGA;
   run = RA;
   detector = DC;

}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // Check if the ion chamber is included in the detector construction.
  if(detector->Get_Ion_Cham_Flag() == "on")
    {
      ion_cham_layer_hc_id = SDman->GetCollectionID("ion_cham_layer/HitsCollection");
      ion_cham_shell_hc_id = SDman->GetCollectionID("ion_cham_shell/HitsCollection");
    }

  // Check if the sample is included in the detector construction.
  if(detector->Get_Sample_Flag() == "on")
     {
	   sample_hc_id = SDman->GetCollectionID("sample/HitsCollection");
     }

  // The experimental hall must be included by definition so no need to test this.
  exp_hall_hc_id = SDman->GetCollectionID("exp_hall/HitsCollection");

  // Increment the event counter by 1.
  event_id = evt->GetEventID() + 1;

  // Print the event number every 10000 lines.
  if (event_id < 101 || event_id%10000 == 0)
  G4cout << "-------------------------------------------->>> Event "
	 << event_id << G4endl;

}

void EventAction::EndOfEventAction(const G4Event* evt) {
	G4HCofThisEvent * HCE = evt->GetHCofThisEvent();

	DetectorHitsCollection* ion_cham_layer_hc = 0; // Initialise a pointer for the actual hits collection for the ion chamber layers.
	DetectorHitsCollection* ion_cham_shell_hc = 0; // Initialise a pointer for the actual hits collection for the ion chamber shell.
	DetectorHitsCollection* exp_hall_hc = 0; // Initialise a pointer for the actual hits collection for the experimental hall.

	// Number of interactions in each sensitive volume.
	G4int num_trigs_ion_cham_layer = 0; // Counter for the number of interactions in the ion chamber layers.
	G4int num_trigs_ion_cham_shell = 0; // Counter for the number of interactions in the ion chamber shell.
	G4int num_trigs_ion_cham_all = 0;   // Counter for the number of interactions in the ion chamber overall.
	G4int num_trigs_exp_hall = 0; // Counter for the number of interactions in the experimental hall.

	// Init. counter for total number of interactions in a detector that can trigger the readout.
	G4int tot_num_trigs = 0;

	// Init. vector to store all hits.
	vector<DetectorHits*> ion_cham_all_hits;

	// Init. vector to store photon hits.
	vector<DetectorHits*> ion_cham_phot_hits;

	// Init. vector to store all hits.
	vector<G4double> all_hit_times;

	// Init. vector to store photon hits.
	vector<G4double> phot_hit_times;

	// Init. vector to store all hits in time ordered vector.
	vector<DetectorHits*> time_sort_all_hits;

	// Init. vector to store photon hits in time ordered vector.
	vector<DetectorHits*> time_sort_phot_hits;

	if (HCE) {

		if (detector->Get_Ion_Cham_Flag() == "on") {

			// Get the hits collection ID for the gas layers of the ion chamber.
			ion_cham_layer_hc = (DetectorHitsCollection*) (HCE->GetHC(ion_cham_layer_hc_id));

            // Test if there are interaction entries to read out.
			if (ion_cham_layer_hc->entries() > 0) {

				// Assume that the layer that of the ionisation chamber gas stack is the middle layer.
				sens_layer = ceil((float)detector->Get_Num_Gas_Layers() / 2.0);

				// Fix the volume that is the sensitive volume.
				sens_vol = "ion_cham_layer_param_phys";

				// Only need to process the event if energy is deposited between the plates - so check this.
                if (Energy_Dep_Between_Plates(ion_cham_layer_hc)){


					// Add the hits to the vector.
					Vectorize_Hits(ion_cham_all_hits, ion_cham_phot_hits, ion_cham_layer_hc, all_hit_times, phot_hit_times);
					num_trigs_ion_cham_layer += 1;
					num_trigs_ion_cham_all += 1;
					tot_num_trigs += 1;


					// Get the hits collection ID for the shell of the ion chamber.
					ion_cham_shell_hc = (DetectorHitsCollection*) (HCE->GetHC(ion_cham_shell_hc_id));


					// Test if there are interaction entries to read out.
					if (ion_cham_shell_hc->entries() > 0){

						// Add the hits to the vector.
						Vectorize_Hits(ion_cham_all_hits, ion_cham_phot_hits, ion_cham_shell_hc, all_hit_times, phot_hit_times);
						num_trigs_ion_cham_shell += 1;
						num_trigs_ion_cham_all += 1;
						tot_num_trigs += 1;

					}
                }
			}
		}
	}

	// Initialise the variables for the current event.
	energy_keV_dep_plates_prim = 0.0;  // Total energy deposited  for primary photon interaction between the ionisation plates for the current event.
	energy_keV_dep_plates_scat = 0.0;  // Total energy deposited  for scattered photon interaction between the ionisation plates for the current event.
	energy_keV_dep_plates_fluor = 0.0; // Total energy deposited  for scattered photon interaction between the ionisation plates for the current event.

	// Test if there have been any interactions in any of the layers of the ion chamber.
    if (num_trigs_ion_cham_layer > 0){

     	// Now time sort the vector of all interactions.
	    Time_Sort(ion_cham_all_hits, all_hit_times, time_sort_all_hits);

	    // Now time sort the vector of photon interactions.
	    Time_Sort(ion_cham_phot_hits, phot_hit_times, time_sort_phot_hits);

	    // Process the interactions in the event.
	    Process_Ion_Cham_Interacs(time_sort_all_hits, time_sort_phot_hits);
	}

}

bool EventAction::Energy_Dep_Between_Plates(DetectorHitsCollection* hits_col)
{
	// Loop over the interactions in the hits collection.
	for (G4int i =0; i< (G4int)hits_col->GetSize(); i++){

		// Get the hit (interaction.)
		DetectorHits* cur_hit = (*hits_col)[i];

		// Test to see if the central ion chamber layer (assumed to be the one with the collection electrodes) has an interaction in it.
		if ((unsigned int)cur_hit->GetCopyNumber() == sens_layer){

			// Check that energy was actually deposited, i.e. it wasn't just a Rayleigh scatter.
			// Get the energyof the interaction.
		    G4double energy_keV = abs(cur_hit->GetEnergyDep() / keV);
		    G4bool is_rayleigh = Is_Rayleigh(energy_keV);
		    if (not is_rayleigh){

		    	// Found a valid energy deposit between the electrodes.
		    	return true;
		    }
		}
	}
	// Haven't managed to find a valid energy deposit.
	return false;
}

void EventAction::Process_Ion_Cham_Interacs(vector<DetectorHits*> &time_sort_all_hits, vector<DetectorHits*>	&time_sort_phot_hits)
{


//   G4cout << "-------------- NEW EVENT -----------------------" << G4endl;
//
//   for (unsigned int i = 0; i < time_sort_all_hits.size(); i++) {
//
//	   // Get the current hit.
//	   DetectorHits* ion_cham_hit = time_sort_all_hits[i];
//
//	   G4cout << ion_cham_hit->GetTrackID() << " "
//	          << ion_cham_hit->GetParentID() << " "
//	          << ion_cham_hit->GetEnergyDep() / keV << " "
//	          << ion_cham_hit->GetGlobalTime() / ns << " "
//	          << ion_cham_hit->GetParticle() << " "
//	          << ion_cham_hit->GetCopyNumber() << " "
//	          << ion_cham_hit->GetPreProcess() << " "
//	          << ion_cham_hit->GetVolume() << " "
//	          << G4endl;
//   }

   // Check if there is a photon interaction stored.
   if ((unsigned int)time_sort_phot_hits.size() > 0){

	   // Get the primary photon interaction.
	   DetectorHits* first_phot_interac = time_sort_phot_hits[0];

       // Test if the photon if a primary or fluorescence photon.
	   G4bool is_fluor  = Is_Fluor(first_phot_interac->GetTrackID(), first_phot_interac->GetParentID());
	   if (not is_fluor){

		   // It is primary!

		   // Test if the interaction happened between the plates.
		   if (((unsigned int)first_phot_interac->GetCopyNumber() == sens_layer)
			   && ((G4String)first_phot_interac->GetVolume() == sens_vol)){

			   // The first interaction happened between the plates.

			   // Now, need to determine if all of the charge is collected at the electrode.

			   // Sum all energies between plates and assign to the total for the primary interaction occurring between the plates.
			   Sum_Energies(time_sort_all_hits, energy_keV_dep_plates_prim, sens_layer, sens_vol);
			   Dump_Photon_Interac_Energy(energy_keV_dep_plates_prim,  1);

		   } else {

			   // The first interaction happened somewhere else not between the plates which must have been a scatter.

			   // Sum all energies between plates and assign to the total for the primary interaction being a scatter.
			   Sum_Energies(time_sort_all_hits, energy_keV_dep_plates_scat, sens_layer, sens_vol);
			   Dump_Photon_Interac_Energy(energy_keV_dep_plates_scat,  2);
		   }
	   }
   }
}

void EventAction::Sum_Energies(vector<DetectorHits*> &all_hits, G4double &total_energy_keV, unsigned int layer, G4String volume){

	// Loop over the interactions and grab hist collection with sorted time stamp.
	for (G4int i = 0; i < (G4int)all_hits.size(); i++){

		// Get the hit (interaction.)
		DetectorHits* cur_hit = all_hits[i];

		// Test if the layer of the interaction is between the plates.
		if (((unsigned int)cur_hit->GetCopyNumber() == layer) && ((G4String)cur_hit->GetVolume() == volume)){

			// Add the energy
			total_energy_keV += abs(cur_hit->GetEnergyDep() / keV);

		}
	}
}

void EventAction::Dump_Photon_Interac_Energy(const G4double dump_energy_keV,  const unsigned int tag){

	// Get the output file pointer.
  	out_file_ptr = run->Get_File_Ptr();

  	// Write the position to file.
  	*out_file_ptr << tag << " " << dump_energy_keV << G4endl;
}

void EventAction::Dump_Prim_Photon_Interac_Pos(const G4double x,  const G4double y, const G4double z){

	// Get the output file pointer.
  	out_file_ptr = run->Get_File_Ptr();

  	// Write the position to file.
  	*out_file_ptr << x << " " << y << " " << z << G4endl;
}

bool EventAction::Is_Photon(const G4String particle_ref){

    // Determine whether or not the photon is a fluorescence or not (primary).
    if (particle_ref == "gamma"){
    	return true;
    }
    else if (particle_ref == "e-"){
    	// It must be an electron as we only deal with photons and electrons in these sims.
    	return false;
    }
    else{
    	// Safety net in case something strange happens
    	abort();
    }
}

bool EventAction::Is_Rayleigh(G4double energy_keV){

	// Determine whether or not the photon is a fluorescence or not (primary).
	if (energy_keV > 0.000000001){
    	return false;
    }
    return true;
}

bool EventAction::Is_Photoelectron(G4String pre_process){

	// Determine whether or not the electron is from the initial primary photon interaction.
	if (pre_process == "phot"){
    	return true;
    }
    return false;
}

bool EventAction::Is_Fluor(const unsigned int track, const unsigned int parent){

	// Determine whether or not the photon is a fluorescence or not (primary).
    if (track == 1 && parent == 0){
    	return false;
    }
    return true;
}

void EventAction::Vectorize_Hits(vector<DetectorHits*>	&ion_cham_all_hits, vector<DetectorHits*> &ion_cham_phot_hits, DetectorHitsCollection* hits_col, vector<G4double> &all_hit_times, vector<G4double> &phot_hit_times)
{
	// Loop over the interactions in the hits collection.
	for (G4int i =0; i< (G4int)hits_col->GetSize(); i++){

		// Get the hit (interaction.)
	    DetectorHits* cur_hit = (*hits_col)[i];

	    // Append the hit to the vector.
	    ion_cham_all_hits.push_back(cur_hit);

		// Check if the interaction is a photon.
		G4bool is_photon = Is_Photon(cur_hit->GetParticle());
		if (is_photon){
			ion_cham_phot_hits.push_back(cur_hit);

			// Append the time stamp to the vector.
			phot_hit_times.push_back((*hits_col)[i]->GetGlobalTime());
		}

		// Append the time stamp to the vector.
		all_hit_times.push_back((*hits_col)[i]->GetGlobalTime());
	}
}


void EventAction::Time_Sort(vector<DetectorHits*>	&all_hits, vector<G4double> &hit_times, vector<DetectorHits*> &time_sort_all_hits)
{
    // Generate a vector of sorted interaction time stamps.
	vector<G4double> sorted_times = hit_times;

    // Sort the vector of times.
	sort(sorted_times.begin(), sorted_times.end());

	// Loop over the interactions and grab hist collection with sorted time stamp.
	for (G4int i = 0; i < (G4int)sorted_times.size(); i++){

		// Now loop again to encompass the un-sorted times.
	 	for(G4int j = 0; j < (G4int)sorted_times.size(); j++){

	 		// Test if time stamp of current hit is the required one.
	 	    if (hit_times.at(j) == sorted_times.at(i)){

	 	    	// Get the hit (interaction.)
	 	    	DetectorHits* sort_hit = all_hits[j];

	 	    	// Append the hit to the vector.
	 	    	time_sort_all_hits.push_back(sort_hit);
	 		    break;
	 	    }
	 	}
	}
}

