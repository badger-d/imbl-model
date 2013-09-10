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
  virtual void Set_Ion_Cham_Flag(G4String val){ion_cham_flag = val;};
  virtual void Set_Sample_Flag(G4String val){sample_flag = val;};
  
  virtual G4String Get_Ion_Cham_Flag(){return ion_cham_flag;};
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
  G4LogicalVolume* ion_cham_log;
  G4VPhysicalVolume* exp_hall_log;
  G4LogicalVolume* sample_log;

  // Physical volume pointers.
  G4VPhysicalVolume* ion_cham_phys;
  G4VPhysicalVolume* exp_hall_phys;
  G4VPhysicalVolume* sample_phys;

  // Sensitive detector pointers.
  SensitiveDet* exp_hall_sd; // Declare pointer to experimental hall.
  SensitiveDet* ion_cham_sd; // Declare pointer to ion chamber sensitive detector.
  SensitiveDet* sample_sd;   // Declare pointer to sample sensitive detector.

  // Material pointers
  G4Material* hall_mat;              // Declare experimental hall material.
  G4Material* ion_cham_fill_gas_mat; // Declare ion chamber fill gas material.
  G4Material* sample_mat;            // Declare sample material.

  // Detector messenger pointer.
  DetectorMessenger* detectorMessenger;

  // Flags to switch components on/off.
  G4String sample_flag;
  G4String ion_cham_flag;
  
  // Variable that will be modified by messenger file.
  G4double sample_x;


};

#endif
