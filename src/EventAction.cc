#include "G4Event.hh"
#include "EventAction.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

using namespace std;

EventAction::EventAction(RunAction* RA, DetectorConstruction* DC)
{
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
	vector<DetectorHits*> all_hits;

	// Init. vector to store photon hits.
	vector<DetectorHits*> phot_hits;

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

				// Assume that the layer that of the ionisation chamber gas stack is the middle layer (numbering starts from 0, so sens_layer is 1 for a 3 layer volume.).
				sens_layer = floor((float)detector->Get_Num_Gas_Layers() / 2.0);

				// Fix the volume that is the sensitive volume.
				sens_vol = "ion_cham_layer_param_phys";

				// Fix the volume that is the shell volume.
				shell_vol = "ion_cham_shell_phys";

				// Only need to process the event if energy is deposited between the plates - so check this.
				if (Energy_Dep_Between_Plates(ion_cham_layer_hc)){

					// Add the hits to the vector.
					Vectorize_Hits(all_hits, phot_hits, ion_cham_layer_hc, all_hit_times, phot_hit_times);
					num_trigs_ion_cham_layer += 1;
					num_trigs_ion_cham_all += 1;
					tot_num_trigs += 1;

					// Get the hits collection ID for the shell of the ion chamber.
					ion_cham_shell_hc = (DetectorHitsCollection*) (HCE->GetHC(ion_cham_shell_hc_id));

					// Test if there are interaction entries to read out.
					if (ion_cham_shell_hc->entries() > 0){

						// Add the hits to the vector.
						Vectorize_Hits(all_hits, phot_hits, ion_cham_shell_hc, all_hit_times, phot_hit_times);
						num_trigs_ion_cham_shell += 1;
						num_trigs_ion_cham_all += 1;
						tot_num_trigs += 1;
					}

					// Get the hits collection ID for the experimental hall.
					exp_hall_hc = (DetectorHitsCollection*) (HCE->GetHC(exp_hall_hc_id));

					// Test if there are interaction entries to read out.
					if (exp_hall_hc->entries() > 0){

						// Add the hits to the vector.
						Vectorize_Hits(all_hits, phot_hits, exp_hall_hc, all_hit_times, phot_hit_times);
						num_trigs_exp_hall += 1;
						tot_num_trigs += 1;
					}
				}
			}
		}
	}

	// Initialise the variables for the current event.
	energy_keV_dep_plates_prim_vola = 0.0;  // Total energy deposited for primary photon interaction between the ionisation plates for the current event.
	energy_keV_dep_plates_prim_volb = 0.0;  // Total energy deposited for primary photon interaction in regions other than between the ionisation plates for the current event.
	energy_keV_dep_plates_scat_vola = 0.0;  // Total energy deposited for scattered photon interaction between the ionisation plates for the current event.
	energy_keV_dep_plates_fluor_vola = 0.0; // Total energy deposited for scattered photon interaction between the ionisation plates for the current event.
	energy_keV_dep_plates_eloss_vola = 0.0; // Total energy deposited for photons whose progeny reached the collecting electrode - volume A.
	energy_keV_dep_plates_eloss_volb = 0.0; // Total energy deposited for photons whose progeny reached the collecting electrode - volume A.

	// Test if there have been any interactions in the sensitive volume of the ion chamber.
	if (num_trigs_ion_cham_layer > 0){

		// Now time sort the vector of all interactions.
		Time_Sort(all_hits, all_hit_times, time_sort_all_hits);

		// Now time sort the vector of photon interactions.
		Time_Sort(phot_hits, phot_hit_times, time_sort_phot_hits);

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
		// G4cout <<  (unsigned int)cur_hit->GetCopyNumber() << " " << sens_layer << G4endl;
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

void EventAction::Process_Ion_Cham_Interacs(vector<DetectorHits*> &time_sort_all_hits, vector<DetectorHits*> &time_sort_phot_hits)
{

	// Init. output flag for file write.
	int output_flag = -1;

	// Init. var for first interaction being a scatter of any kind.
	bool has_scat = false;

	// Init. var for a photon interaction being a primary or fluorescence.
	bool is_primary = false;

	// Init. var for fluorescence photons being present.
	bool has_fluor = false;

	// Check if there is a photon interaction stored.
	if ((unsigned int)time_sort_phot_hits.size() > 0){

		// There is a photon in the list.

		// Get the primary photon interaction.
		DetectorHits* first_phot_hit = time_sort_phot_hits[0];

		// Test if the photon if a primary or fluorescence photon.
		is_primary = Is_Primary(first_phot_hit->GetTrackID(), first_phot_hit->GetParentID());

		// Init. vector to store the progeny of the first interaction.
		vector<DetectorHits*> first_phot_progeny;

		// Init. vector to store the progeny of fluorescence interacitons.
		vector<DetectorHits*> fluor_phot_progeny;

		if (is_primary){

			// The photon is a primary.

			// A strange type of interaction where the photon is not recorded may happen, so avert this.
			unsigned int interaction_index = Get_First_Photon_Index(time_sort_all_hits);

			// Add this first interaction to the progeny for the energy sum.
			first_phot_progeny.push_back(time_sort_all_hits[interaction_index]);

			// Now check if the photon has scattered.
			has_scat = Has_Scat(time_sort_phot_hits);

			// Some of the energy must have left the sensitive volume otherwise its indistinguishable from the events for E_{ap}
			//bool has_left_sens = Has_Left_Sens(time_sort_phot_hits, sens_layer, sens_vol);

			if (has_scat){

				// -------------- Process for $k_{asc}$. --------------------- //

				// Sum the energies for the progeny.
				// Just search through all interactions, do not group according to tracks.
				energy_keV_dep_plates_scat_vola = Sum_Energies_In_Hit_Vol(time_sort_all_hits, sens_layer, sens_vol);

				// Dump the sum energy to file.
				output_flag = 2;
				Dump_Photon_Interac_Energy(energy_keV_dep_plates_scat_vola, output_flag);

			} else{

				// It is a primary photon interaction that has not scattered.

				// Get the photon progeny knowing that the interaction was a photoelectric absorption.
				Get_Photon_Progeny_Phot(time_sort_all_hits, interaction_index, first_phot_progeny);

				// Sum the energies for the progeny
				energy_keV_dep_plates_prim_vola = Sum_Energies_In_Hit_Vol(first_phot_progeny, sens_layer, sens_vol);

				// Dump the sum energy to file.
				output_flag = 0;
				Dump_Photon_Interac_Energy(energy_keV_dep_plates_prim_vola, output_flag);

				// Sum the energies for the progeny
				energy_keV_dep_plates_prim_volb = Sum_Energies_Not_In_Hit_Vol(first_phot_progeny, sens_layer, sens_vol, shell_vol);

				// Dump the sum energy to file.
				output_flag = 1;
				Dump_Photon_Interac_Energy(energy_keV_dep_plates_prim_volb, output_flag);

				/*
				G4cout << energy_keV_dep_plates_prim_vola << " " << energy_keV_dep_plates_prim_volb << G4endl;

				bool non_layer = false;

				// Loop over the interactions to find tracks with the corresponding parent ID.
				for (unsigned int i = 0; i < (unsigned int)time_sort_all_hits.size(); i++){

					// Get the hit (interaction).
					DetectorHits* cur_hit = time_sort_all_hits[i];

					if (cur_hit->GetVolume() != "ion_cham_layer_param_phys"){

						non_layer = true;

					}

				}

				if (non_layer){

					G4cout << "-------------- NEW EVENT --------------" << G4endl;
					G4cout << energy_keV_dep_plates_prim_vola << " " << energy_keV_dep_plates_prim_volb << G4endl;
					for (unsigned int j = 0; j < (unsigned int)time_sort_all_hits.size(); j++){

						// Get the hit (interaction).
						DetectorHits* cur_hit = time_sort_all_hits[j];

						G4cout << cur_hit->GetTrackID() << " " <<  cur_hit->GetParentID() << " " << abs(cur_hit->GetEnergyDep()) << " " << cur_hit->GetProcess() << " " << cur_hit->GetParticle() << " " <<  cur_hit->GetVolume() << " " << (unsigned int)cur_hit->GetCopyNumber() <<G4endl;
					}

				}
				 */


				// Now search for fluorescence photons.
				has_fluor = Has_Fluor(time_sort_phot_hits);

				//				if ((energy_keV_dep_plates_prim_vola > 0.0000001) && (energy_keV_dep_plates_prim_volb > 0.0000001)){
				//					G4cout << time_sort_phot_hits.size() << " " << has_fluor<< " " << energy_keV_dep_plates_prim_vola + energy_keV_dep_plates_prim_volb << " " << energy_keV_dep_plates_prim_vola << " " << energy_keV_dep_plates_prim_volb << G4endl;
				//				}

				if (has_fluor){

					// Dig out the progeny of the fluorescence photons.

					// Init. vector to store the indices of the fluorescence interactions.
					vector<unsigned int> fluor_phot_indices;

					// Get the indices of the fluorescence photons.
					Get_Fluor_Photon_Indices(time_sort_all_hits, fluor_phot_indices);

					// Get the fluorescence photon progeny.
					Get_Photon_Progeny_Fluor(time_sort_all_hits, fluor_phot_indices, fluor_phot_progeny);

					// Sum the energies for the progeny
					energy_keV_dep_plates_fluor_vola = Sum_Energies_In_Hit_Vol(fluor_phot_progeny, sens_layer, sens_vol);

					// Dump the sum energy to file.
					output_flag = 3;
					Dump_Photon_Interac_Energy(energy_keV_dep_plates_fluor_vola, output_flag);

				}
			}
		}
	}
}


bool EventAction::Has_Fluor(vector<DetectorHits*> &phot_hits){

	// Determine if a fluorescence photon is present.
	for (int i = 0; i < (int)phot_hits.size(); i++){

		// Get the hit (interaction).
		DetectorHits* cur_hit = phot_hits[i];

		// Test if the photon is a primary or fluorescence photon.
		bool is_primary  = EventAction::Is_Primary(cur_hit->GetTrackID(), cur_hit->GetParentID());
		if (not is_primary){

			// A fluorescence photon was found.
			return true;
		}
	}

	return false;
}


bool EventAction::Has_Scat(vector<DetectorHits*> &phot_hits){

	if ((int)phot_hits.size() <= 1){

		// It is not a scatter.
		return false;
	}

	// It might be a scatter but it might be a photo plus fluor, so filter.

	// Init. a var. for counting the number of primary photon interactions.
	G4int prim_phot_interacs = 0;

	// Determine if the photon has scattered.
	for (int i = 0; i < (int)phot_hits.size(); i++){

		// Get the hit (interaction).
		DetectorHits* cur_hit = phot_hits[i];

		// Test if the photon is a primary or fluorescence photon.
		bool is_primary  = EventAction::Is_Primary(cur_hit->GetTrackID(), cur_hit->GetParentID());
		if (is_primary){

			// Increment the counter by one.
			prim_phot_interacs += 1;
		}

	}

	if (prim_phot_interacs >= 2){

		// It is a scatter.
		return true;
	}

	return false;
}

bool EventAction::Has_Left_Sens(vector<DetectorHits*> &phot_hits, unsigned int layer, G4String volume){

	// Determine if the photon has scattered.
	for (int i = 0; i < (int)phot_hits.size(); i++){

		// Get the hit (interaction).
		DetectorHits* cur_hit = phot_hits[i];

		// Test if the layer of the interaction is between the plates.
		if (not (((unsigned int)cur_hit->GetCopyNumber() == layer) && ((G4String)cur_hit->GetVolume() == volume))){

			return true;
		}
	}
	return false;
}


bool EventAction::Has_Relevant_Parent(unsigned int parent, vector<unsigned int> parent_vector){

	for (int i = 0; i < (int)parent_vector.size(); i++){

		if (parent == (unsigned int)parent_vector[i]) {

			return true;
		}
	}
	return false;
}

void EventAction::Add_Track_ID(unsigned int track, vector<unsigned int> &parent_vector){

	// Assume that the track ID is not in the master parent ID list.
	bool in_parent_vector = false;

	for (int i = 0; i < (int)parent_vector.size(); i++){

		if (track == parent_vector[i]){

			in_parent_vector = true;
			break;
		}

	}

	if (not in_parent_vector) {

		parent_vector.push_back(track);
	}
}

unsigned int EventAction::Get_First_Photon_Index(vector<DetectorHits*> &all_hits){

	// Get the hit index that corresponds to the first photon interaction.

	// Assume there are no photons.
	bool is_photon = false;

	// Init. index for photon interaction.
	unsigned int i;

	// Loop over the interactions.
	for (i = 0; i < (unsigned int)all_hits.size(); i++){

		// Get the interaction that is the start point.
		DetectorHits* cur_hit = all_hits[i];

		// Test to see if the interaction is a photon.
		is_photon = EventAction::Is_Photon(cur_hit->GetParticle());
		if (is_photon){

			break;
		}
	}

	// A photon must be found.
	assert(is_photon == true);

	return i;
}

void EventAction::Get_Photon_Progeny_Phot(vector<DetectorHits*> &all_hits, unsigned int init_interac_index, vector<DetectorHits*> &daughter_tracks){

	// Layer store check.
	vector<unsigned int> parent_vector;

	// Get the interaction that is the start point.
	DetectorHits* first_interac = all_hits[init_interac_index];

	// The first interaction must be a photon.
	assert(first_interac->GetParticle() == "gamma");

	// Init. var for maximum interaction to loop up to.
	unsigned int max = (unsigned int)all_hits.size();

	// Just a normal photoelectric interaction, so will loop over all interactions.
	max = (unsigned int)all_hits.size();

	// Store the track ID.
	parent_vector.push_back(first_interac->GetTrackID());

	// Loop over the interactions to find tracks with the corresponding parent ID.
	for (unsigned int i = init_interac_index + 1; i < max; i++){

		// Get the hit (interaction).
		DetectorHits* cur_hit = all_hits[i];

		// Get the current parent.
		unsigned int cur_parent = cur_hit->GetParentID();

		// Determine if the interaction has a relevant parent.
		bool has_relevant_parent = Has_Relevant_Parent(cur_parent, parent_vector);

		if (has_relevant_parent){

			// Add the hit to the vector.
			daughter_tracks.push_back(cur_hit);

			// Get the current track.
			unsigned int cur_track = cur_hit->GetTrackID();

			// Add this track ID as a valid parent ID for the mother interaction.
			Add_Track_ID(cur_track, parent_vector);
		}
	}
}

void EventAction::Get_Photon_Progeny_Scat(vector<DetectorHits*> &all_hits, unsigned int init_interac_index, vector<DetectorHits*> &daughter_tracks){

	// Layer store check.
	vector<unsigned int> parent_vector;

	// Get the interaction that is the start point.
	DetectorHits* first_interac = all_hits[init_interac_index];

	// The first interaction must be a photon.
	assert(first_interac->GetParticle() == "gamma");

	// Init. var for maximum interaction to loop up to.
	unsigned int max = (unsigned int)all_hits.size();

	// Init. loop var.
	unsigned int i;

	// Shorten the stack to the interaction before the scatter.
	for (i = init_interac_index + 1; i < (unsigned int)all_hits.size(); i++){

		// Get the hit (interaction).
		DetectorHits* cur_hit = all_hits[i];

		// Determine if the interaction is a primary.
		bool is_primary = Is_Primary(cur_hit->GetTrackID(), cur_hit->GetParentID());

		if (is_primary){
			break;

		}

	}

	// Set the maximum interaction to loop up to to be the break point.
	max = i;

	// Store the track ID.
	parent_vector.push_back(first_interac->GetTrackID());

	// Loop over the interactions to find tracks with the corresponding parent ID.
	for (unsigned int i = init_interac_index + 1; i < max; i++){

		// Get the hit (interaction).
		DetectorHits* cur_hit = all_hits[i];

		// Get the current parent.
		unsigned int cur_parent = cur_hit->GetParentID();

		// Determine if the interaction has a relevant parent.
		bool has_relevant_parent = Has_Relevant_Parent(cur_parent, parent_vector);

		if (has_relevant_parent){

			// Add the hit to the vector.
			daughter_tracks.push_back(cur_hit);

			// Get the current track.
			unsigned int cur_track = cur_hit->GetTrackID();

			// Add this track ID as a valid parent ID for the mother interaction.
			Add_Track_ID(cur_track, parent_vector);
		}
	}
}

void EventAction::Get_Fluor_Photon_Indices(vector<DetectorHits*> &all_hits, vector<unsigned int> &indices){

	// Get the hit indices that corresponds to fluorescence photons.

	// Loop over the interactions.
	for (unsigned int i = 0; i < (unsigned int)all_hits.size(); i++){

		// Get the interaction that is the start point.
		DetectorHits* cur_hit = all_hits[i];

		// Test to see if the interaction is a photon.
		bool is_photon = EventAction::Is_Photon(cur_hit->GetParticle());

		if (is_photon){

			// Test to see if the photon is a primary or fluorescence photon.
			bool is_primary = EventAction::Is_Primary(cur_hit->GetTrackID(), cur_hit->GetParentID());

			if (not is_primary) {

				// It is a fluorescence photon so add the index to the vector.
				indices.push_back(i);
			}
		}
	}
}

void EventAction::Get_Photon_Progeny_Fluor(vector<DetectorHits*> &all_hits, vector<unsigned int> &indices, vector<DetectorHits*> &daughter_tracks){

	// Init. interaction index.
	unsigned int fluor_interac_index;

	// Layer store check.
	vector<unsigned int> parent_vector;

	// Loop over the fluorescence interaciton indices.
	for (unsigned int i = 0; i < (unsigned int)indices.size(); i++){

		// Pass the current fluorescence index.
		fluor_interac_index = indices[i];

		// Get the hit (interaction).
		DetectorHits* fluor_photon_hit = all_hits[fluor_interac_index];

		// Add the parent hit to the vector.
		daughter_tracks.push_back(fluor_photon_hit);

		// Store the track ID.
		parent_vector.push_back(fluor_photon_hit->GetTrackID());

		// Add one to the interaction index to skip on to progeny.
		fluor_interac_index += 1;

		// Loop over the interactions to find the fluorescence progeny.
		for (unsigned int j = fluor_interac_index; j < (unsigned int)all_hits.size(); j++){

			// Get the hit (interaction).
			DetectorHits* cur_hit = all_hits[j];

			// Get the current parent.
			unsigned int cur_parent = cur_hit->GetParentID();

			// Determine if the interaction has a relevant parent.
			bool has_relevant_parent = Has_Relevant_Parent(cur_parent, parent_vector);

			if (has_relevant_parent){

				// Add the hit to the vector.
				daughter_tracks.push_back(cur_hit);

				// Get the current track.
				unsigned int cur_track = cur_hit->GetTrackID();

				// Add this track ID as a valid parent ID for the mother interaction.
				Add_Track_ID(cur_track, parent_vector);

			}
		}
	}
}

bool EventAction::All_Electron_Interacs_In_Sens(vector<DetectorHits*> &all_hits){

	// Loop through the interactions.
	for (G4int i = 0; i < (G4int)all_hits.size(); i++){

		// Get the hit (interaction.)
		DetectorHits* cur_hit = all_hits[i];

		// Get the process that caused the interaction.
		G4String pre_process = cur_hit->GetPreProcess();

		// Find if the interaction is an electron.
		bool is_photoelectron = EventAction::Is_Photoelectron(pre_process);

		// Test if the interaction is an electron.
		if (is_photoelectron){

			// Test if the interaction did not stay between the plates.
			if (((unsigned int)cur_hit->GetCopyNumber() != sens_layer)
					|| ((G4String)cur_hit->GetVolume() != sens_vol)){
				return false;
			}
		}

	}
	return true;
}

G4double EventAction::Sum_Energies_In_Hit_Vol(vector<DetectorHits*> &all_hits, unsigned int sensitive_layer, G4String sensitive_volume){

	// Sum the energies of the interactions that ARE INSIDE the designate volume and layer.
	G4double total_energy_keV = 0.0;

	// Loop over the interactions and grab hist collection with sorted time stamp.
	for (G4int i = 0; i < (G4int)all_hits.size(); i++){

		// Get the hit (interaction.)
		DetectorHits* cur_hit = all_hits[i];

		// Test if the layer of the interaction is between the plates.
		if (((unsigned int)cur_hit->GetCopyNumber() == sensitive_layer) && ((G4String)cur_hit->GetVolume() == sensitive_volume)){

			// Add the energy
			total_energy_keV += abs(cur_hit->GetEnergyDep() / keV);

		}
	}
	return total_energy_keV;
}

G4double EventAction::Sum_Energies_Not_In_Hit_Vol(vector<DetectorHits*> &all_hits, unsigned int sensitive_layer, G4String sensitive_volume, G4String shell_volume){

	// Sum the energies of the interactions that ARE NOT INSIDE the designate volume and layer.
	G4double total_energy_keV = 0.0;

	// Loop over the interactions and grab hist collection with sorted time stamp.
	for (G4int i = 0; i < (G4int)all_hits.size(); i++){

		// Get the hit (interaction.)
		DetectorHits* cur_hit = all_hits[i];

		// Test if the layer of the interaction is between the plates.
		if ((G4String)cur_hit->GetVolume() == sensitive_volume){

			if ((unsigned int)cur_hit->GetCopyNumber() != sensitive_layer){

				total_energy_keV += abs(cur_hit->GetEnergyDep() / keV);

			}
		}
		else{

			// It has to be in the shell.
			if ((G4String)cur_hit->GetVolume() == shell_volume){
				total_energy_keV += abs(cur_hit->GetEnergyDep() / keV);

			}
		}
	}
	return total_energy_keV;
}


void EventAction::Dump_Photon_Interac_Energy(const G4double dump_energy_keV,  const unsigned int output_flag){

	if (dump_energy_keV > 0.000000001){

		// Get the output file pointer.
		out_file_ptr = run->Get_File_Ptr();

		// Write the position to file.
		*out_file_ptr << output_flag << " " << dump_energy_keV << G4endl;

	}
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

bool EventAction::Is_Primary(const unsigned int track, const unsigned int parent){

	// Determine whether or not the photon is a primary or not (fluorescence).
	if (track == 1 && parent == 0){
		return true;
	}
	return false;
}

void EventAction::Vectorize_Hits(vector<DetectorHits*>	&all_hits, vector<DetectorHits*> &phot_hits, DetectorHitsCollection* hits_col, vector<G4double> &all_hit_times, vector<G4double> &phot_hit_times)
{
	// Loop over the interactions in the hits collection.
	for (G4int i =0; i< (G4int)hits_col->GetSize(); i++){

		// Get the hit (interaction.)
		DetectorHits* cur_hit = (*hits_col)[i];

		// Append the hit to the vector.
		all_hits.push_back(cur_hit);

		// Append the time stamp to the vector.
		all_hit_times.push_back((*hits_col)[i]->GetGlobalTime());

		// Check if the interaction is a photon.
		G4bool is_photon = Is_Photon(cur_hit->GetParticle());
		if (is_photon){
			phot_hits.push_back(cur_hit);

			// Append the time stamp to the vector.
			phot_hit_times.push_back((*hits_col)[i]->GetGlobalTime());
		}
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

