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
	vector<DetectorHits*>	ion_cham_hits;

	// Init. vector to store all hits.
	vector<G4double> hit_times;

	// Init. vector to store all hits in time ordered vector.
	vector<DetectorHits*>	time_sort_all_hits;

	if (HCE) {

		if (detector->Get_Ion_Cham_Flag() == "on") {

			// Get the hits collection ID for the gas layers of the ion chamber.
			ion_cham_layer_hc = (DetectorHitsCollection*) (HCE->GetHC(ion_cham_layer_hc_id));

            // Test if there are interaction entries to read out.
			if (ion_cham_layer_hc->entries() > 0) {

				// Add the hits to the vector.
                Vectorize_Hits(ion_cham_hits, ion_cham_layer_hc, hit_times);
                num_trigs_ion_cham_layer += 1;
                num_trigs_ion_cham_all += 1;
                tot_num_trigs += 1;
			}

			// Get the hits collection ID for the shell of the ion chamber.
  		    ion_cham_shell_hc = (DetectorHitsCollection*) (HCE->GetHC(ion_cham_shell_hc_id));

  		    // Test if there are interaction entries to read out.
  		  	if (ion_cham_shell_hc->entries() > 0){

          	    // Add the hits to the vector.
                Vectorize_Hits(ion_cham_hits, ion_cham_shell_hc, hit_times);
                num_trigs_ion_cham_shell += 1;
                num_trigs_ion_cham_all += 1;
                tot_num_trigs += 1;

			}
		}
	}

	// Test if there have been any interactions in any part of the ion chamber.
    if (num_trigs_ion_cham_all > 0){

     	// Now time sort the vector of all interactions.
	    Time_Sort(ion_cham_hits, hit_times, time_sort_all_hits);

	    // Process the interactions in the event.
	    Process_Ion_Cham_Interacs(ion_cham_hits)

//	    for (G4int i = 0; i < time_sort_all_hits.size(); i++) {
//	        DetectorHits* sort_hit = time_sort_all_hits[i];
//	        sort_hit->GetDeltaEnergy();
//	        G4cout << sort_hit->GetDeltaEnergy() << G4endl;
//	    }
	}
}

void EventAction::Process_Ion_Cham_Interacs(vector<DetectorHits*>	&all_hits)
{
    G4cout << "------------- HIT ION CHAMBER -------------" << G4endl;
}


void EventAction::Vectorize_Hits(vector<DetectorHits*>	&all_hits, DetectorHitsCollection* hits_col, vector<G4double>	&hit_times)
{
	// Loop over the interactions in the hits collection.
	for (G4int i =0; i< (G4int)hits_col->GetSize(); i++){

		// Get the hit (interaction.)
	    DetectorHits* cur_hit = (*hits_col)[i];

	    // Append the hit to the vector.
		all_hits.push_back(cur_hit);

		// Append the time stamp to the vector.
		hit_times.push_back((*hits_col)[i]->GetGlobalTime());
	}
}


//void EventAction::Time_Sort(vector<DetectorHits*>	&all_hits, vector<G4double> &hit_times, vector<DetectorHits*>* time_sort_all_hits)
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

