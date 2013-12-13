//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file field/field03/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id$
// 

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include <vector>

class CalorimeterSD;
class FieldSetup;

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class SensitiveDet;
class ExpPhantomSD;
class DetectorInfo;
class G4Box;

class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
    ~DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);          
      
     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);
     
     virtual G4VPhysicalVolume* Construct();

     void Update_Geometry();
     
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

     // Set the ion chamber properties based on the name selected by the user.
     void Set_Ion_Cham_Properties();
  
     void PrintParameters();

  private:
     
     G4UniformMagField*    fMagField;         // pointer to the magnetic field
     FieldSetup*        emFieldSetup;

     void DefineMaterials();
     void ComputeCalorParameters();
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
      G4Material* ion_cham_sens_mat;     // Define the material for the sensitive volume of the ion chamber.

      // Detector messenger pointer.
      DetectorMessenger* detector_messenger;

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
      G4double ion_cham_wall_dz_mm;
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
      std::vector<G4double> ion_cham_gas_dz_mm;         // Store for the thicknesses of the gas layers along z-axis.

      // Specify the correction value that ensures there are no boundary clashes.
      G4double correc_fac;


};

#endif
