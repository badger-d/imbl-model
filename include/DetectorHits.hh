//DetectorHits.hh
//-------------------------------------


#ifndef DetectorHits_h
#define DetectorHits_h

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class DetectorHits : public G4VHit
{
public:

  DetectorHits();
  virtual ~DetectorHits();

  virtual void Print();

  virtual void SetParticle(G4String val){particle = val;};
  virtual void SetTrackID(G4int val){trackID = val;};
  virtual void SetParentID(G4int val){parentID = val;};
  virtual void SetPosition(G4ThreeVector val){position = val;};
  virtual void SetDeltaPosition(G4ThreeVector val){delta_position = val;};
  virtual void SetDeltaEnergy(G4double val){delta_energy = val;};
  virtual void SetEnergyDep(G4double val){energy_dep = val;};
  virtual void SetEnergyKinetic(G4double val){energy_kinetic = val;};
  virtual void SetDeltaMomentum(G4ThreeVector val){delta_momentum = val;};
  virtual void SetVolume(G4String val){volume = val;};
  virtual void SetProcess(G4String val){process = val;};
  virtual void SetGlobalTime(G4double val){time = val;}
  virtual void SetGlobalTranslation(G4ThreeVector val){global_trans = val;};
  virtual void SetLocalTranslation(G4ThreeVector val){local_trans = val;};
  virtual void SetCopyNumber(G4int val){copy_no = val;};
  virtual void SetMultiplicity(G4int val){multiplicity = val;};
  virtual void SetStepLength(G4double val){step_length = val;};
  virtual void SetPreProcess(G4String val){preprocess = val;};
  virtual void SetOrigLog(G4String val){orig_log = val;};
  virtual void SetOrigPos(G4ThreeVector val){orig_pos = val;};

  virtual G4int GetMultiplicity(){return multiplicity;};
  virtual G4int GetCopyNumber(){return copy_no;};
  virtual G4ThreeVector GetGlobalTranslation(){return global_trans;};
  virtual G4ThreeVector GetLocalTranslation(){return local_trans;};
  virtual G4String GetVolume(){return volume;};
  virtual G4String GetParticle(){return particle;};
  virtual G4int GetTrackID(){return trackID;};
  virtual G4int GetParentID(){return parentID;};
  virtual G4ThreeVector GetPosition(){return position;};
  virtual G4ThreeVector GetDeltaPosition(){return delta_position;};
  virtual G4double GetDeltaEnergy(){return delta_energy;};
  virtual G4double GetEnergyDep(){return energy_dep;};
  virtual G4double GetEnergyKinetic(){return energy_kinetic;};
  virtual G4ThreeVector GetDeltaMomentum(){return delta_momentum;};
  virtual G4String GetProcess(){return process;};
  virtual G4double GetGlobalTime(){return time;};
  virtual G4double GetStepLength(){return step_length;};
  virtual G4String GetPreProcess(){return preprocess;};
  virtual G4String GetOrigLog(){return orig_log;};
  virtual G4ThreeVector GetOrigPos(){return orig_pos;};

private:
  G4int multiplicity;
  G4int copy_no;
  G4String orig_log;
  G4ThreeVector orig_pos;
  G4ThreeVector global_trans;
  G4ThreeVector local_trans;
  G4double time;
  G4String volume;
  G4ThreeVector delta_momentum;
  G4double delta_energy;
  G4double energy_dep;
  G4double energy_kinetic;
  G4ThreeVector delta_position;
  G4ThreeVector position;
  G4int parentID;
  G4int trackID;
  G4String particle;
  G4String process;
  G4String preprocess;
  G4double step_length;
};

//----------------------------------------------------------------------------
typedef G4THitsCollection<DetectorHits> DetectorHitsCollection;

extern G4Allocator<DetectorHits> DetectorHitsAllocator;

//----------------------------------------------------------------------------
#endif


