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

  if(detector->Get_Detector_Flag() == "on")
    {
      detectorHC_ID = SDman->GetCollectionID("detector/HitsCollection");
    }

  if(detector->Get_Sample_Flag() == "on")
     {
	   sampleHC_ID = SDman->GetCollectionID("sample/HitsCollection");
     }

  hallHC_ID = SDman->GetCollectionID("hall/HitsCollection");


  G4int eventID = evt->GetEventID() + 1;

  if (eventID < 101 || eventID%10000 == 0)
  G4cout << "-------------------------------------------->>> Event "
	 << eventID << G4endl;

}

void EventAction::EndOfEventAction(const G4Event* evt) {
	G4HCofThisEvent * HCE = evt->GetHCofThisEvent();

	DetectorHitsCollection* detector_HC = 0;
	DetectorHitsCollection* hall_HC = 0;

	// Number of interactions in each sensitive volume.
	G4int numTrigs_detector = 0;
	G4int numTrigs_hall = 0;

	// Total number of interactions.
	G4int numTrigs = 0;

	G4int eventID = evt->GetEventID() + 1;

	if (HCE) {

		if (detector->Get_Detector_Flag() == "on") {
			detector_HC = (DetectorHitsCollection*) (HCE->GetHC(detectorHC_ID));
			if (detector_HC->entries() > 0) {
				numTrigs_detector += detector_HC->entries();
				numTrigs++;
			}
		}

		hall_HC = (DetectorHitsCollection*) (HCE->GetHC(hallHC_ID));

		if (hall_HC->entries() > 0) {
			numTrigs_hall += hall_HC->entries();
		}


	}

    if (numTrigs > 0) {

		// Get the output file pointer.
		outFile_ptr = run->Get_File_Ptr();

		// Output the detector data.
		if (numTrigs_detector > 0) {
			for (G4int i = 0; i < numTrigs_detector; i++) {
				G4int partrefnum = 0;
				G4String partref = (*detector_HC)[i]->GetParticle();

				if (partref == "electron") {
					partrefnum = 0;
				} else if (partref == "gamma") {
					partrefnum = 1;
				}

				*outFile_ptr << (*detector_HC)[i]->GetTrackID() << " "
						<< (*detector_HC)[i]->GetParentID() << " "
						<< (*detector_HC)[i]->GetEnergyDep() / keV << " "
						<< -(*detector_HC)[i]->GetDeltaEnergy() / keV << " "
						<< (*detector_HC)[i]->GetPosition().x() / cm << " "
						<< (*detector_HC)[i]->GetPosition().y() / cm << " "
						<< (*detector_HC)[i]->GetPosition().z() / cm << " "
						<< (*detector_HC)[i]->GetGlobalTime() / ns << " "
						<< (*detector_HC)[i]->GetPreProcess() << " " \
						<< partrefnum << " "
						<< i << " "
						<< 1000 << " "
						<< (*detector_HC)[i]->GetProcess() << " "
						<< eventID << " "
						<< G4endl;
			}
		}

	}
}

