#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SensitiveDet.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4ios.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4BREPSolidPolyhedra.hh"
#include "G4RotationMatrix.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include <vector>
#include <assert.h>


DetectorConstruction::DetectorConstruction()
{
  detectorMessenger = new DetectorMessenger(this);

  // Select the ion chamber to be included / excluded.
  ion_cham_flag = "on";
  Set_Ion_Cham_Flag(ion_cham_flag);

  // Select sample to be included / excluded.
  sample_flag = "off";
  Set_Sample_Flag(sample_flag);

  // Sensitive volume pointers.
  exp_hall_sd = 0; // Initialize pointer to experimental hall sensitive detector.
  ion_cham_sd = 0; // Initialize pointer to ion chamber sensitive detector.
  
  // Initialize the thickness of the sample.
  sample_x = 2.0 * mm;

}

DetectorConstruction::~DetectorConstruction()
{
   delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
   DefineMaterials();
   return ConstructGeometry();
}

// Materials definitions

void DetectorConstruction::DefineMaterials()
{
   G4double a;                     // Material atomic mass.
   G4double z;                     // Material atomic number.
   G4double density;               // Material density.
   G4String name;                  // Material name.
   G4String symbol;                // Elemental symbol.
   G4int ncomponents;              // Number of constituents.
   G4double fractionmass;          // Fractional mass.
   G4double temperature, pressure; // Temperature and pressure.

   // Define elements that might be required in the model.

   G4Element* elH  = new G4Element("Hydrogen", "H", 1. ,  1.01*g/mole);

   G4Element* elHe = new G4Element("Helium", "He", 2., 4.002*g/mole);

   G4Element* elC = new G4Element("Carbon", "C", 6., 12.01*g/mole);

   G4Element* elN = new G4Element("Nitrogen", "N", 7., 14.01*g/mole);

   G4Element* elO = new G4Element("Oxygen", "O", 8., 16.00*g/mole);

   G4Element* elF = new G4Element("Flourine", "F", 9., 18.99*g/mole);

   G4Element* elNe = new G4Element("Neon", "Ne", 10., 20.1797*g/mole);

   G4Element* elAl = new G4Element("Aluminium","Al", 13., 26.98*g/mole);

   G4Element* elSi = new G4Element("Silicon", "Si", 14., 28.09*g/mole);

   G4Element* elAr =  new G4Element("Argon", "Ar", 18., 39.948*g/mole);

   G4Element* elCa = new G4Element("Calcium", "Ca", 20., 40.08*g/mole);

   G4Element* elCr = new G4Element("Chromium", "Cr", 24., 52.00*g/mole);

   G4Element* elFe = new G4Element("Iron", "Fe", 26., 55.847*g/mole);

   G4Element* elNi = new G4Element("Nickel","Ni", 28., 58.70*g/mole);

   G4Element* elCu = new G4Element("Copper", "Cu", 29., 63.546*g/mole);

   G4Element* elKr =  new G4Element("Krypton", "Kr", 36., 83.798*g/mole);

   G4Element* elXe = new G4Element("Xenon", "Xe", 54., 131.293*g/mole);

   G4Element* elPb = new G4Element("Lead", "Pb", 82., 207.19*g/mole);

   // The materials are defined below to allow impurities to be added.

   // Air
   G4Material* Air = new G4Material(name="Air", 0.001290*g/cm3, ncomponents=8);
   Air->AddElement(elC, fractionmass=0.0124*perCent);
   Air->AddElement(elN, fractionmass=75.87697130*perCent);
   Air->AddElement(elO, fractionmass=23.1781*perCent);
   Air->AddElement(elNe, fractionmass=0.0018*perCent);
   Air->AddElement(elHe, fractionmass=0.00052*perCent);
   Air->AddElement(elAr, fractionmass=0.93*perCent);
   Air->AddElement(elXe, fractionmass=0.0000087*perCent);
   Air->AddElement(elKr, fractionmass=0.0002*perCent);

   // Vacuum:
   density = universe_mean_density;
   pressure = 3.e-18*pascal;
   temperature = 273*kelvin;
   new G4Material(name="Vacuum", z=1., a=1.01*g/mole, density, kStateGas, pressure, temperature);

   // Specify the default materials used in the experimental hall.
   ion_cham_fill_gas_mat = Air;
   hall_mat = Air;
   sample_mat = Air;
}


// Geometry definitions

G4VPhysicalVolume* DetectorConstruction::ConstructGeometry()
{
    G4cout << "Constructing the geometry..." << G4endl;

    // Experimental hall (world volume)

    G4double exp_hall_x = .6*m; // World depth.
    G4double exp_hall_y = .5*m; // World width.
    G4double exp_hall_z = .5*m; // World height.

    G4Box* exp_hall_box
      = new G4Box("exp_hall_box", exp_hall_x * 0.5, exp_hall_y * 0.5, exp_hall_z * 0.5);

    G4LogicalVolume* exp_hall_log
      = new G4LogicalVolume(exp_hall_box, hall_mat,"nexp_hall_log",0,0,0);

    G4VPhysicalVolume* exp_hall_phys
      = new G4PVPlacement(0, G4ThreeVector(), "nexp_hall", exp_hall_log, 0, false, 0);

    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // Make experimental hall a sensitive volume

    if (! exp_hall_sd)
      {
    	exp_hall_sd = new SensitiveDet("/detector/hall");
        SDman->AddNewDetector(exp_hall_sd);
      }
    exp_hall_log->SetSensitiveDetector(exp_hall_sd);

    exp_hall_log->SetVisAttributes (G4VisAttributes::Invisible);


    // Switch the detector on/off

    if(Get_Ion_Cham_Flag() == "on")
      {

        G4cout << "Including the ion chamber .........................................................." << G4endl;

        // Detector
        G4double detector_x = 1.0*mm;
        G4double detector_y = 10.0*mm;
        G4double detector_z = 10.0*mm;

        G4ThreeVector ion_cham_pos = G4ThreeVector(20.0, 0.0, 0.0);

        G4Box* detector = new G4Box("nDetector", detector_x * 0.5, detector_y * 0.5, detector_z * 0.5);
        ion_cham_log = new G4LogicalVolume(detector, ion_cham_fill_gas_mat,"nDetector_log",0,0,0);

        ion_cham_phys = new G4PVPlacement(0,                // Rotation.
        								  ion_cham_pos,     // (x,y,z).
        								  ion_cham_log,     // Logical volume.
        								  "ion_cham_phys", // Name.
        								  exp_hall_log,      // Mother volume.
        								  false,            // No boolean operations.
        								  0);               // Copy number.


        G4VisAttributes* detector_visAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));

        // Make PD a sensitive volume
        ion_cham_log->SetVisAttributes(detector_visAtt);

        if (! ion_cham_sd)
          {
        	ion_cham_sd = new SensitiveDet("/detectors/detector");
            SDman->AddNewDetector(ion_cham_sd);
          }
        ion_cham_log->SetSensitiveDetector(ion_cham_sd);

      }
    else
      {
        G4cout << "Excluding the ion chamber .........................................................." << G4endl;
      }


    if(Get_Sample_Flag() == "on")
      {

        G4cout << "Including the sample .................................................." << G4endl;


		G4double sample_y = 10.0 * mm;
		G4double sample_z = 10.0 * mm;

		G4ThreeVector sample_pos = G4ThreeVector(10.0, 0.0, 0.0);

		G4VSolid* sample = new G4Box("nSample", Get_Sample_X() * 0.5, sample_y * 0.5, sample_z * 0.5);

		sample_log = new G4LogicalVolume(sample, sample_mat, "sample_log");

		sample_phys = new G4PVPlacement(0,              // Rotation.
										sample_pos,     // (x,y,z).
										sample_log,     // Logical volume.
										"sample_phys", // Name.
										exp_hall_log,    // Mother volume.
										false,          // No boolean operations.
										0);             // Copy number.


        if(! sample_sd)
			{
        	sample_sd = new SensitiveDet("/detectors/sample");
              SDman->AddNewDetector(sample_sd);
            }

        sample_log->SetSensitiveDetector(sample_sd);

        // Set visualization attributes.
        G4VisAttributes* sample_visAtt = new G4VisAttributes(G4Colour(1.0,0.7,1.0));
        sample_log->SetVisAttributes(sample_visAtt);

      }
    else
    {
    	G4cout << "Excluding the sample .........................................................." << G4endl;
    }



  G4cout << "Detector construction complete" << G4endl;

  return exp_hall_phys;
}

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructGeometry());
}



