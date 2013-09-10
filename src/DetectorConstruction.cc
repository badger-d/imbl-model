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

  // Select detector
  detector_flag = "on";
  Set_Detector_Flag(detector_flag);

  // Select sample
  sample_flag = "on";
  Set_Sample_Flag(sample_flag);

  // Seneitive volume pointers.
  aHallSD = 0;
  aDetectorSD = 0;
  aSampleSD = 0;
  
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
   G4int ncomponents, natoms;      // Number of constituents and atoms.
   G4double fractionmass;          // Fractional mass.
   G4double temperature, pressure; // Temperature and pressure.

   // Elements

   a = 1.01*g/mole;
   G4Element* elH = new G4Element(name="Hydrogen",symbol="H", z =1., a);

   a = 12.01*g/mole;
   G4Element* elC = new G4Element(name="Carbon",symbol="C", z =6., a);

   a = 14.01*g/mole;
   G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z = 7., a);

   a = 16.00*g/mole;
   G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);

   a = 18.998*g/mole;
   G4Element* elF = new G4Element(name="Flourine", symbol="F", z=9., a);

   a = 26.98*g/mole;
   G4Element* elAl = new G4Element(name="Aluminium",symbol="Al", z =13., a);

   a = 28.09*g/mole;
   G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);

   a = 40.08*g/mole;
   G4Element* elCa = new G4Element(name="Calcium", symbol="Ca", z=20., a);

   a = 52.00*g/mole;
   G4Element* elCr = new G4Element(name="Chromium", symbol="Cr", z=24., a);

   a = 55.847*g/mole;
   G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z=26., a);

   a = 58.70*g/mole;
   G4Element* elNi = new G4Element(name="Nickel",symbol="Ni",z=28., a);

   a = 63.546*g/mole;
   G4Element* elCu = new G4Element(name="Copper", symbol="Cu", z=29. ,a);

   a = 207.19*g/mole;
   G4Element* elPb = new G4Element(name="Lead", symbol="Pb", z=82., a);

   // Single element materials

   density = 7.15 *g/cm3;
   G4Material* Cr = new G4Material(name="Cr", density, ncomponents=1);
   Cr->AddElement(elCr, fractionmass=1.);

   density =  7.874*g/cm3;
   G4Material* Fe = new G4Material(name="Fe", density, ncomponents=1);
   Fe->AddElement(elFe, fractionmass=1.);

   density = 8.912 *g/cm3;
   G4Material* Ni = new G4Material(name="Ni", density, ncomponents=1);
   Ni->AddElement(elNi, fractionmass=1.);

   density = 8.920 *g/cm3;
   G4Material* Cu = new G4Material(name="Cu", density, ncomponents=1);
   Cu->AddElement(elCu, fractionmass=1.);

   density =  2.33*g/cm3;
   G4Material* Si = new G4Material(name="Si", density, ncomponents=1);
   Si->AddElement(elSi, fractionmass=1.);

   density = 1.55*g/cm3;
   G4Material* Ca = new G4Material(name="Ca", density, ncomponents=1);
   Ca->AddElement(elCa, fractionmass=1.);

   density = 2.7*g/cm3;
   G4Material* Al = new G4Material(name="Al", density, ncomponents=1);
   Al->AddElement(elAl, natoms=1);

   density = 11.35*g/cm3;
   G4Material* Pb = new G4Material(name="Pb", density, ncomponents=1);
   Pb->AddElement(elPb, fractionmass=1.);

   // Air
   density = .001290*g/cm3;
   G4Material* Air = new G4Material(name="Air", density, ncomponents=2);
   Air->AddElement(elN, fractionmass=70*perCent);
   Air->AddElement(elO, fractionmass=30*perCent);

   // B-100 Bone equivalent
   density = 1.450*g/cm3;
   G4Material* B100
	= new G4Material(name= "B-100 Bone equivalent", density, ncomponents=6);
   B100->AddElement(elH, fractionmass=.065473);
   B100->AddElement(elC, fractionmass=.536942);
   B100->AddElement(elN, fractionmass=.021500);
   B100->AddElement(elO, fractionmass=.032084);
   B100->AddElement(elF, fractionmass=.167415);
   B100->AddElement(elCa, fractionmass=.176585);
   B100->GetIonisation()->SetMeanExcitationEnergy(85.9*eV);

   // Vacuum:
   density = universe_mean_density;
   pressure = 3.e-18*pascal;
   temperature = 273*kelvin;

   new G4Material(name="Vacuum", z=1., a=1.01*g/mole, density, kStateGas, pressure, temperature);

   // Specify the default materials.
   detector_mat = Si;
   hall_mat = Air;
   sample_mat = B100;

}


// Geometry definitions

G4VPhysicalVolume* DetectorConstruction::ConstructGeometry()
{
    G4cout << "Constructing the geometry..." << G4endl;

    // Experimental hall (world volume)

    G4double expHall_x = .6*m; // World depth.
    G4double expHall_y = .5*m; // World width.
    G4double expHall_z = .5*m; // World height.

    G4Box* expHall_box
      = new G4Box("expHall_box", expHall_x * 0.5, expHall_y * 0.5, expHall_z * 0.5);

    G4LogicalVolume* expHall_log
      = new G4LogicalVolume(expHall_box, hall_mat,"nExpHall_log",0,0,0);

    G4VPhysicalVolume* expHall_phys
      = new G4PVPlacement(0, G4ThreeVector(), "nExphall", expHall_log, 0, false, 0);

    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // Make experimental hall a sensitive volume

    if (! aHallSD)
      {
        aHallSD = new SensitiveDet("/detector/hall");
        SDman->AddNewDetector(aHallSD);
      }
    expHall_log->SetSensitiveDetector(aHallSD);

    expHall_log->SetVisAttributes (G4VisAttributes::Invisible);


    // Switch the detector on/off

    if(Get_Detector_Flag() == "on")
      {

        G4cout << "Including the detector .........................................................." << G4endl;

        // Detector
        G4double detector_x = 1.0*mm;
        G4double detector_y = 10.0*mm;
        G4double detector_z = 10.0*mm;

        G4ThreeVector detector_pos = G4ThreeVector(20.0, 0.0, 0.0);

        G4Box* detector = new G4Box("nDetector", detector_x * 0.5, detector_y * 0.5, detector_z * 0.5);
        detector_log = new G4LogicalVolume(detector, detector_mat,"nDetector_log",0,0,0);

        detector_phys = new G4PVPlacement(0,                // Rotation.
        								  detector_pos,     // (x,y,z).
        								  detector_log,     // Logical volume.
        								  "nDetector_phys", // Name.
        								  expHall_log,      // Mother volume.
        								  false,            // No boolean operations.
        								  0);               // Copy number.


        G4VisAttributes* detector_visAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));

        // Make PD a sensitive volume
        detector_log->SetVisAttributes(detector_visAtt);

        if (! aDetectorSD)
          {
            aDetectorSD = new SensitiveDet("/detectors/detector");
            SDman->AddNewDetector(aDetectorSD);
          }
        detector_log->SetSensitiveDetector(aDetectorSD);

      }
    else
      {
        G4cout << "Excluding the detector .........................................................." << G4endl;
      }


    if(Get_Sample_Flag() == "on")
      {

        G4cout << "Including the sample .................................................." << G4endl;


		G4double sample_y = 10.0 * mm;
		G4double sample_z = 10.0 * mm;

		G4ThreeVector sample_pos = G4ThreeVector(10.0, 0.0, 0.0);

		G4VSolid* sample = new G4Box("nSample", Get_Sample_X() * 0.5, sample_y * 0.5, sample_z * 0.5);

		sample_log = new G4LogicalVolume(sample, sample_mat, "nSample_log");

		sample_phys = new G4PVPlacement(0,              // Rotation.
										sample_pos,     // (x,y,z).
										sample_log,     // Logical volume.
										"nSample_phys", // Name.
										expHall_log,    // Mother volume.
										false,          // No boolean operations.
										0);             // Copy number.


        if(! aSampleSD)
			{
        	  aSampleSD = new SensitiveDet("/detectors/sample");
              SDman->AddNewDetector(aSampleSD);
            }

        sample_log->SetSensitiveDetector(aSampleSD);

        // Set visualization attributes.
        G4VisAttributes* sample_visAtt = new G4VisAttributes(G4Colour(1.0,0.7,1.0));
        sample_log->SetVisAttributes(sample_visAtt);

      }
    else
    {
    	G4cout << "Excluding the sample .........................................................." << G4endl;
    }



  G4cout << "Detector construction complete" << G4endl;

  return expHall_phys;
}

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructGeometry());
}



