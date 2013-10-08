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
class ElectricFieldSetup;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include <vector>
#include "G4ElectricField.hh"
#include "G4FieldManager.hh"

using namespace std;

class DetectorConstruction : public G4VUserDetectorConstruction
{

public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  // Set inclusion of the beam-line components [on, off]
  virtual void Set_Ion_Cham_Flag(G4String val){ion_cham_flag = val;};
  virtual void Set_Sample_Flag(G4String val){sample_flag = val;};
  
  // Get inclusion of the beam-line components [on, off]
  virtual G4String Get_Ion_Cham_Flag(){return ion_cham_flag;};
  virtual G4String Get_Sample_Flag(){return sample_flag;};

  // Set the sample position.
  virtual void Set_Sample_Pos(G4ThreeVector val){sample_pos = val;}; // Collimator frame depth.

  // Get the sample position.
  virtual G4ThreeVector Get_Sample_Pos(){return sample_pos;};

  // Set the number of gas layers in the ionisation chamber.
  virtual void Set_Num_Gas_Layers(unsigned int val){num_gas_layers = val;}; // Collimator frame depth.

  // Get the number of gas layers in the ionisation chamber.
  virtual unsigned int Get_Num_Gas_Layers(){return num_gas_layers;};

  // Control updating of the geometry.
  virtual void Update_Geometry(); // Updates geometry.

  // Set the ion chamber properties based on the name selected by the user.
  void Set_Ion_Cham_Properties();

private:
  void DefineMaterials();
  G4VPhysicalVolume* Construct_Geometry();

  // Logical volume pointers.
  G4LogicalVolume* ion_cham_shell_log;     // Declare pointer to logical volume of shell of ion chamber.
  G4LogicalVolume* ion_cham_sens_log;      // Declare pointer to logical volume of sensitive volume of ion chamber.
  G4LogicalVolume* ion_cham_layer_log;     // Declare pointer to logical volume of arbitrary layer of ion chamber volume.
  G4LogicalVolume* exp_hall_log;
  G4LogicalVolume* sample_log;


  // Physical volume pointers.
  G4VPhysicalVolume* ion_cham_shell_phys;       // Declare pointer to physical volume of shell of ion chamber.
  G4VPhysicalVolume* ion_cham_sens_phys;        // Declare pointer to physical volume of sensitive volume of ion chamber.
  G4VPhysicalVolume* ion_cham_layer_param_phys; // Declare pointer to physical volume of stack parameterisation of ion chamber gases.
  G4VPhysicalVolume* exp_hall_phys;
  G4VPhysicalVolume* sample_phys;

  // Sensitive detector pointers.
  SensitiveDet* exp_hall_sd;       // Declare pointer to experimental hall.
  SensitiveDet* ion_cham_layer_sd; // Declare pointer to ion chamber layer sensitive detectors.
  SensitiveDet* ion_cham_shell_sd; // Declare pointer to ion chamber shell sensitive detector.
  SensitiveDet* sample_sd;         // Declare pointer to sample sensitive detector.

  // Material pointers
  G4Material* exp_hall_mat;          // Declare experimental hall material.
  G4Material* ion_cham_shell_mat;    // Define the material for the shell of the ion chamber.
  G4Material* sample_mat;            // Declare sample material.
  G4Material* vacuum_mat;            // Declare vacuum material.

  // Detector messenger pointer.
  DetectorMessenger* detector_messenger;

  // Electromagnetic field pointer.
  ElectricFieldSetup* em_field_setup;

  // Flags to switch components on/off.
  G4String sample_flag;
  G4String ion_cham_flag;
  
  // Variables that will be modified by messenger file.
  G4ThreeVector ion_cham_pos; // Position of the ion chamber.
  G4ThreeVector sample_pos;   // Position of the sample.
  G4String ion_cham_name;     // Name of the required ion chambe.

  // Specify the dimensions of the external section of the Al box that is the shell of the ion chamber.
  G4ThreeVector ion_cham_ext_mm;

  // Specify the dimensions of the internal section of the Al box that is the shell of the ion chamber.
  // The shell thickness is 10.0*mm all the way around.
  G4double ion_cham_wall_thick_mm;
  G4ThreeVector ion_cham_int_mm;

  // Specify the dimensions of the sensitive volume of the ion chamber.
  G4ThreeVector ion_cham_sens_mm;

  // Specify the dimensions of the entrance / exit windows.
  G4double ion_cham_ent_win_rad_mm;  // Radius of entrance window in mm.
  G4double ion_cham_exit_win_rad_mm; // Radius of exit window in mm.

  // Generic cylinder parameters.
  G4double gener_cyl_r;  // Inner radius..
  G4double gener_cyl_SP; // Tube segment.
  G4double gener_cyl_EP; // Delta angle.

  // Properties of ion chamber gas layers.
  std::vector<G4String> ion_cham_gas_mat_name_store;   // Store for the fill-gas name.
  std::vector<G4ThreeVector> ion_cham_gas_shape_store; // Store for the dimensions of the fill gas layer.
  std::vector<G4ThreeVector> ion_cham_gas_pos_store;   // Store for the locations of the fill gas layers.
  std::vector<G4Material*> ion_cham_gas_mat_store;     // Store for the pointer to the fill gas material object.
  unsigned int num_gas_layers;                         // Number of layers of gas.
  std::vector<G4double> ion_cham_gas_thick_mm;         // Store for the thicknesses of the gas layers along z-axis.

  // Specify the correction value that ensures there are no boundary clashes.
  G4double correc_fac;

  // Electric field for volume of the ionisation chamber.
  G4ElectricField* ic_em_field;
  G4FieldManager* local_ic_field_mgr;
};

#endif
