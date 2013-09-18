#include "SensitiveDet.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4VProcess.hh"

SensitiveDet::SensitiveDet(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="HitsCollection");
  HCID = -1;
}
SensitiveDet::~SensitiveDet()
{
}

void SensitiveDet::Initialize(G4HCofThisEvent* HCE)
{
  hitsCollection = new DetectorHitsCollection
                          (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
   }

  HCE->AddHitsCollection(HCID,hitsCollection);
}

G4bool SensitiveDet::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

   particle = aStep->GetTrack()->GetDefinition()->GetParticleName();
   parentID = aStep->GetTrack()->GetParentID();
   G4ThreeVector delta_energy(aStep->GetDeltaEnergy());

   G4String process = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
   if(process == "Transportation") {

	   return false;
   }

   DetectorHits* newHit = new DetectorHits();
   newHit->SetGlobalTime(aStep->GetPostStepPoint()->GetGlobalTime());
   newHit->SetVolume(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());
   newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
   newHit->SetParentID(parentID);
   newHit->SetParticle(particle);
   newHit->SetStepLength(aStep->GetTrack()->GetStepLength());
   newHit->SetPosition(aStep->GetTrack()->GetPosition());
   newHit->SetDeltaPosition(aStep->GetDeltaPosition());
   newHit->SetDeltaEnergy(aStep->GetDeltaEnergy());
   newHit->SetEnergyDep(aStep->GetTotalEnergyDeposit());
   newHit->SetEnergyKinetic(aStep->GetPostStepPoint()->GetKineticEnergy());
   newHit->SetDeltaMomentum(delta_energy);
   newHit->SetProcess(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
   newHit->SetGlobalTranslation(aStep->GetPreStepPoint()->GetTouchable()->GetTranslation());
   newHit->SetLocalTranslation(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation());
   newHit->SetCopyNumber(aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
   newHit->SetPreProcess("None");

     const G4VProcess* proc = aStep->GetTrack()->GetCreatorProcess();
     if (proc != 0)
     {
    	 const G4String pre_process = proc->GetProcessName();
		 newHit->SetPreProcess(pre_process);
     }

     const G4LogicalVolume* logVol = aStep->GetTrack()->GetLogicalVolumeAtVertex();
     const G4String orig_log = logVol->GetName();
     newHit->SetOrigLog(orig_log);

     const G4ThreeVector& orig_pos = aStep->GetTrack()->GetVertexPosition();
     newHit->SetOrigPos(orig_pos);

   hitsCollection->insert(newHit);

   return true;

}

void SensitiveDet::EndOfEvent(G4HCofThisEvent*)
{
    if (verboseLevel>0)
    {
      G4int NbHits = hitsCollection->entries();
        G4cout << "\n-------->Hits Collection: in this event there are " << NbHits
               << " hits in " << GetVolume() << ": " << G4endl;
        for (G4int i=0;i<NbHits;i++) (*hitsCollection)[i]->Print();
    }
}

