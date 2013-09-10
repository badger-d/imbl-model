#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <fstream>
#include <vector>

class DetectorConstruction;
class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class HeaderInfo;

using std::ifstream;
using std::ofstream;


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*);
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  virtual void Set_Isotropy(G4String val) {isotropy = val;}
  virtual void Set_Energy(G4double val) {energy = val;}
  virtual void Set_Origin(G4ThreeVector val){origin = val;};

  virtual G4String Get_Isotropy(){return isotropy;};
  virtual G4double Get_Energy(){return energy;};
  virtual G4ThreeVector Get_Origin(){return origin;};

private:

  DetectorConstruction* detector;
  PrimaryGeneratorMessenger* messenger;
  G4ParticleGun* particleGun;

  G4double energy;
  G4String isotropy;
  G4ThreeVector origin;

};

#endif

