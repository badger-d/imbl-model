#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class SensitiveDet;
class ExpPhantomSD;
class DetectorInfo;
class G4Box;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include <vector>

using namespace std;

class DetectorConstruction : public G4VUserDetectorConstruction
{

public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  // Control inclusion of housing structures [on, off]
  virtual void Set_Detector_Flag(G4String val){detector_flag = val;};
  virtual void Set_Sample_Flag(G4String val){sample_flag = val;};
  
  virtual G4String Get_Detector_Flag(){return detector_flag;};
  virtual G4String Get_Sample_Flag(){return sample_flag;};

  // Control the sample thickness.
  virtual void Set_Sample_X(G4double val){sample_x = val;}; // Collimator frame depth.
  virtual G4double Get_Sample_X(){return sample_x;};

  // Control updating of the geometry.
  virtual void UpdateGeometry(); // Updates geometry.

private:
  void DefineMaterials();
  G4VPhysicalVolume* ConstructGeometry();

  // Logical volume pointers.
  G4LogicalVolume* detector_log;
  G4LogicalVolume* sample_log;
  G4VPhysicalVolume* expHall_log;

  // Physical volume pointers.
  G4VPhysicalVolume* detector_phys;
  G4VPhysicalVolume* sample_phys;
  G4VPhysicalVolume* expHall_phys;

  // Sensitive detector pointers.
  SensitiveDet* aDetectorSD;
  SensitiveDet* aSampleSD;
  SensitiveDet* aHallSD;

  // Material pointers
  G4Material* hall_mat;   // Experimental hall material.
  G4Material* detector_mat;   // Collimator material.
  G4Material* sample_mat;   // Collimator backstop material.

  // Detector messenger pointer.
  DetectorMessenger* detectorMessenger;

  // Flags to switch components on/off.
  G4String sample_flag;
  G4String detector_flag;
  
  // Variable that will be modified by messenger file.
  G4double sample_x;


};

#endif
