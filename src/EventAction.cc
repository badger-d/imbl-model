#include "G4Event.hh"
#include "EventAction.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "DetectorConstruction.hh"

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
      ion_cham_hc_id = SDman->GetCollectionID("ion_cham/HitsCollection");
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

	DetectorHitsCollection* ion_cham_hc = 0; // Initialise a pointer for the actual hits collection for the ion chamber.
	DetectorHitsCollection* exp_hall_hc = 0; // Initialise a pointer for the actual hits collection for the experimental hall.

	// Number of interactions in each sensitive volume.
	G4int num_trigs_ion_cham = 0; // Counter for the number of interactions in the ion chamber.
	G4int num_trigs_exp_hall = 0; // Counter for the number of interactions in the experimental hall.

	// Init. counter for total number of interactions in a detector that can trigger the readout.
	G4int tot_num_trigs = 0;

	if (HCE) {

		if (detector->Get_Ion_Cham_Flag() == "on") {

			ion_cham_hc = (DetectorHitsCollection*) (HCE->GetHC(ion_cham_hc_id));

			if (ion_cham_hc->entries() > 0) {
				// If there are triggers in the ion chamber, add them to the count.
				num_trigs_ion_cham += ion_cham_hc->entries();
				tot_num_trigs++;
			}
		}

		exp_hall_hc = (DetectorHitsCollection*) (HCE->GetHC(exp_hall_hc_id));

		if (exp_hall_hc->entries() > 0) {
			num_trigs_exp_hall += exp_hall_hc->entries();
		}


	}

	// If any valid detector has interactions in it the read out the transport.
    if (tot_num_trigs > 0) {

		// Get the output file pointer.
		out_file_ptr = run->Get_File_Ptr();

		// Output the detector data.
		if (num_trigs_ion_cham > 0) {
			for (G4int i = 0; i < num_trigs_ion_cham; i++) {
				G4int partrefnum = 0;
				G4String partref = (*ion_cham_hc)[i]->GetParticle();

				if (partref == "electron") {
					partrefnum = 0;
				} else if (partref == "gamma") {
					partrefnum = 1;
				}

				*out_file_ptr << (*ion_cham_hc)[i]->GetTrackID() << " "
					          << (*ion_cham_hc)[i]->GetParentID() << " "
						      << (*ion_cham_hc)[i]->GetEnergyDep() / keV << " "
						      << -(*ion_cham_hc)[i]->GetDeltaEnergy() / keV << " "
						      << (*ion_cham_hc)[i]->GetPosition().x() / cm << " "
						      << (*ion_cham_hc)[i]->GetPosition().y() / cm << " "
						      << (*ion_cham_hc)[i]->GetPosition().z() / cm << " "
						      << (*ion_cham_hc)[i]->GetGlobalTime() / ns << " "
						      << (*ion_cham_hc)[i]->GetPreProcess() << " "
						      << partrefnum << " "
						      << i << " "
						      << 1000 << " "
						      << (*ion_cham_hc)[i]->GetProcess() << " "
						      << event_id << " "
						      << G4endl;
			}
		}

	}
}

