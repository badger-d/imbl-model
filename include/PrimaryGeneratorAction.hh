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

  // Load the energy array from the file.
  virtual void Load_Energy_Array_Vector();

  // Generate the primary particles.
  virtual void GeneratePrimaries(G4Event*);

  // Set the source emission directional distribution.
  virtual void Set_Isotropy(G4String val) {isotropy = val;}

  // Set the source emission energy.
  virtual void Set_Energy(G4double val) {energy = val;}

  // Set the energy type.
  virtual void Set_Energy_Type(G4String val) {entype = val;}

  // Set the source emission location.
  virtual void Set_Origin(G4ThreeVector val){origin = val;};

  // Set the source diameter.
  virtual void Set_Diameter(G4double val){diameter = val;};

  // Set the source emission position distribution.
  virtual void Set_Distribution(G4String val) {distribution = val;}

  // Get the source emission directional distribution.
  virtual G4String Get_Isotropy(){return isotropy;};

  // Set the source emission energy.
  virtual G4double Get_Energy(){return energy;};

  // Get the source emission location.
  virtual G4ThreeVector Get_Origin(){return origin;};

  // Get the source diameter.
  virtual G4double Get_Diameter(){return diameter;};

  // Get the source emission position distribution.
  virtual G4String Get_Distribution(){return distribution;};

  // Return the type of energy distribution.
  virtual G4String Get_Energy_Type(){return entype;};

private:

  DetectorConstruction* detector;
  PrimaryGeneratorMessenger* messenger;
  G4ParticleGun* particleGun;

  // Get the energy from the array.
  G4double Energy();
  std::vector<double> energyArray;
  int energyOffset;
  G4String entype;
  std::vector<double> probArray;

  // Get the emission location randomly sampled according to the input distribution string.
  G4ThreeVector Distribution(G4String shape);

  // Source particle energy.
  G4double energy;

  // Source emission distribution.
  G4String isotropy;

  // Source emission origin.
  G4ThreeVector origin;

  // Source diameter.
  G4double diameter;

  // String for specifying the emission distribution.
  G4String distribution;

};

#endif

