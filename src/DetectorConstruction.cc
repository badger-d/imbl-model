#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SensitiveDet.hh"
#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "PrimaryGeneratorAction.hh"
#include "StackParameterisation.hh"
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
#include "ElectricFieldSetup.hh"
#include <vector>
#include <assert.h>
#include <map>
#include <numeric>

DetectorConstruction::DetectorConstruction()
{

  // Create an instance of the detector messenger.
  detector_messenger = new DetectorMessenger(this);

  // Create an instance of the electromagnetic field setup.
  em_field_setup = new ElectricFieldSetup();

  // Select the ion chamber to be included / excluded.
  ion_cham_flag = "on";
  Set_Ion_Cham_Flag(ion_cham_flag);

  // Select sample to be included / excluded.
  sample_flag = "off";
  Set_Sample_Flag(sample_flag);

  // Sensitive volume pointers.
  exp_hall_sd = 0;       // Initialise pointer to experimental hall sensitive detector.
  ion_cham_layer_sd = 0; // Initialise pointer to ion chamber layer sensitive detectors.
  ion_cham_shell_sd = 0; // Initialise pointer to ion chamber shell sensitive detector.
  sample_sd = 0;         // Initialise pointer to sample sensitive detector.

  // Initialise the position of the ion chamber.
  ion_cham_pos = G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm);
  
  // Initialise the thickness of the sample.
  sample_pos = G4ThreeVector(2.0*mm, 0.0*mm, 0.0*mm);

  // Specify the correction value that ensures there are no boundary clashes.
  // Typically this is set to 0.01*mm.
  correc_fac = 0.01*mm;

}

DetectorConstruction::~DetectorConstruction()
{
   delete detector_messenger;
   if (em_field_setup) delete em_field_setup;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
   DefineMaterials();
   return Construct_Geometry();
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

   // C
   G4Material* C = new G4Material(name="C", 2.26*g/cm3, ncomponents=1);
   C->AddElement(elC, fractionmass=1.0);

   // Al
   G4Material* Al = new G4Material(name="Al", 2.7*g/cm3, ncomponents=1);
   Al->AddElement(elAl, fractionmass=1.0);

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
   G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a=1.01*g/mole, density, kStateGas, pressure, temperature);

   // Specify the default materials used in the experimental hall.
   ion_cham_shell_mat = Al; // Define the material for the shell of the ion chamber.
   exp_hall_mat = Air;      // Define experimental hall material.
   sample_mat = Air;        // Define sample material.
   vacuum_mat = Vacuum;        // Define vacuum material.
}


// Geometry definitions

G4VPhysicalVolume* DetectorConstruction::Construct_Geometry()
{
    G4cout << "Constructing the geometry..." << G4endl;

    // ------------------------------------- //
    // Experimental hall (world volume).
    // ------------------------------------- //

    G4double exp_hall_x_m = 2.*m; // World depth.
    G4double exp_hall_y_m = 2.*m; // World width.
    G4double exp_hall_z_m = 2.*m; // World height.

    // Define the parameters of the box that is the experimental hall.
    G4Box* exp_hall_box = new G4Box("exp_hall_box", exp_hall_x_m * 0.5, exp_hall_y_m * 0.5, exp_hall_z_m * 0.5);

    // Define the logical volume of the experimental hall.
    exp_hall_log = new G4LogicalVolume(exp_hall_box, exp_hall_mat,"exp_hall_log",0,0,0);

    // Define the physical volume of the experimental hall.
    exp_hall_phys = new G4PVPlacement(0, G4ThreeVector(), "exp_hall", exp_hall_log, 0, false, 0);

    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // Make experimental hall a sensitive volume.

    if (! exp_hall_sd)
      {
    	exp_hall_sd = new SensitiveDet("/detector/exp_hall");
        SDman->AddNewDetector(exp_hall_sd);
      }
    exp_hall_log->SetSensitiveDetector(exp_hall_sd);

    // Set the visualisation attributes for the experimental hall.
    exp_hall_log->SetVisAttributes (G4VisAttributes::Invisible);

    // Generic cylinder parameters used throughout the code.
    gener_cyl_r = 0.0;        // Inner radius..
    gener_cyl_SP = 0.0*deg;   // Tube segment.
    gener_cyl_EP = 360.0*deg; // Delta angle.

    // Switch the ion chamber on/off

    if(Get_Ion_Cham_Flag() == "on")
      {

        G4cout << "Including the ion chamber .........................................................." << G4endl;

        // ------------------------------------- //
        // Ion chamber.
        // ------------------------------------- //

        // Before the stack parameterisation can work, the offsets must be calculated.
        Set_Ion_Cham_Properties();

        // Construct the external box of the ion chamber.
        G4VSolid* ion_cham_shell = new G4Box("ion_cham_shell", ion_cham_ext_mm.getX() * 0.5, ion_cham_ext_mm.getY() * 0.5, ion_cham_ext_mm.getZ() * 0.5);

        // Construct the internal box of the ion chamber.
        G4VSolid* ion_cham_int = new G4Box("ion_cham_int", ion_cham_int_mm.getX() * 0.5, ion_cham_int_mm.getY() * 0.5, ion_cham_int_mm.getZ() * 0.5);

        // Construct the Al shell of the ion chamber.  This will be the subtraction of the internal box from the external box.
        ion_cham_shell = new G4SubtractionSolid("ion_cham_shell", ion_cham_shell, ion_cham_int, 0, G4ThreeVector(0.0*mm, 0.0*mm, 0.0*mm));

        // Cut holes in the front and back walls for the beam to pass through.
        // Define cut in entrance window.
        G4Tubs* ion_cham_ent_cut = new G4Tubs("depdet_be_cut", gener_cyl_r, ion_cham_ent_win_rad_mm, (ion_cham_wall_thick_mm * 0.5) + correc_fac, gener_cyl_SP, gener_cyl_EP);

        // Define cut in entrance window.
        G4Tubs* ion_cham_exit_cut = new G4Tubs("depdet_be_cut", gener_cyl_r, ion_cham_exit_win_rad_mm, (ion_cham_wall_thick_mm * 0.5) + correc_fac, gener_cyl_SP, gener_cyl_EP);

        // Remove the entrance window from the shell.
        ion_cham_shell = new G4SubtractionSolid("ion_cham_shell", ion_cham_shell, ion_cham_ent_cut, 0, G4ThreeVector(0.0*mm, 0.0*mm, -(0.5 * ion_cham_ext_mm.getZ()) + (0.5 * ion_cham_wall_thick_mm)));

        // Remove the exit window from the shell.
        ion_cham_shell = new G4SubtractionSolid("ion_cham_shell", ion_cham_shell, ion_cham_exit_cut, 0, G4ThreeVector(0.0*mm, 0.0*mm, +(0.5 * ion_cham_ext_mm.getZ()) - (0.5 * ion_cham_wall_thick_mm)));

        // Construct the logical volume for the shell of the ion chamber.
        ion_cham_shell_log = new G4LogicalVolume(ion_cham_shell, ion_cham_shell_mat,"ion_cham_shell_log",0,0,0);

        // Construct the physical volume for the shell of the ion chamber.
        ion_cham_shell_phys = new G4PVPlacement(0,                  // Rotation.
        								        ion_cham_pos,       // (x,y,z).
        								        ion_cham_shell_log, // Logical volume.
        								        "ion_cham_phys",    // Name.
        								        exp_hall_log,       // Mother volume.
        								        false,              // No boolean operations.
        								         0);                // Copy number.

        // Now perform the parameterisation of the internal volume of the ion chamber.
        // First, construct a vacuum container that will act as a mother volume for the parameterised layers.
		G4Box* ion_cham_sens = new G4Box("ion_cham_sens", ion_cham_sens_mm.getX() * 0.5,
				                                          ion_cham_sens_mm.getY() * 0.5,
				                                          ion_cham_sens_mm.getZ() * 0.5);

		// Construct the logical volume for the vacuum container that will act as a mother volume for the parameterised layers.
	    ion_cham_sens_log = new G4LogicalVolume(ion_cham_sens, vacuum_mat,"ion_cham_sens_log",0,0,0);

	    // Construct the physical volume for the vacuum container that will act as a mother volume for the parameterised layers.
	    ion_cham_sens_phys = new G4PVPlacement(0,                            // Rotation.
						                       G4ThreeVector(0.0, 0.0, 0.0), // (x,y,z).
						                       ion_cham_sens_log,            // Logical volume.
				        					   "ion_cham__sens_phys",        // Name.
				        					   exp_hall_log,                 // Mother volume.
		            						   false,                        // No boolean operations.
				        					   0);                           // Copy number.

	    // Construct an arbitrary single layer of gas.
	    // The dimensions will be modified by the G4PVParameterised class.
	    G4Box* ion_cham_layer = new G4Box("ion_cham_layer", ion_cham_sens_mm.getX() * 0.5,
	   				                                        ion_cham_sens_mm.getY() * 0.5,
	   				                                        ion_cham_sens_mm.getZ() * 0.5);

	    // Construct the logical volume for the arbitrary single layer of gas.
	    ion_cham_layer_log = new G4LogicalVolume(ion_cham_layer, vacuum_mat, "ion_cham_layer_log",0,0,0);


		// Create each instance of a section of the sensitive volume.

		G4VPVParameterisation* ion_cham_layer_param = new StackParameterisation(num_gas_layers,            // Number of layers.
				                                                                ion_cham_gas_mat_store,    // Vector of material pointers.
				                                                                ion_cham_gas_pos_store,    // Vector of position vectors.
				                                                                ion_cham_gas_shape_store); // Vector of shape vectors.

		// Construct the physical volume for the stack of layers of gas in the ion chamber.
		ion_cham_layer_param_phys = new G4PVParameterised("ion_cham_layer_param_phys", // Name.
				                                          ion_cham_layer_log,          // Logical volume.
		                                                  ion_cham_sens_log,           // Mother logical volume.
		                                                  kZAxis,                      // Axis.
		                                                  num_gas_layers,              // Number of layers.
		                                                  ion_cham_layer_param);       // Parametrisation.

        // Set the visualisation attributes of the ion chamber.
        G4VisAttributes* ion_cham_shell_vis_att = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
        ion_cham_shell_log->SetVisAttributes(ion_cham_shell_vis_att);

        G4VisAttributes* ion_cham_layer_vis_att = new G4VisAttributes(G4Colour(0.3,0.0,1.0));
        ion_cham_layer_log->SetVisAttributes(ion_cham_layer_vis_att);

        // Make the ion chamber layers sensitive volumes.
        if (! ion_cham_layer_sd)
          {
        	ion_cham_layer_sd = new SensitiveDet("/detector/ion_cham_layer");
            SDman->AddNewDetector(ion_cham_layer_sd);
          }

        // Set the logical volume to contain an electric field.
        ion_cham_layer_log->SetFieldManager(em_field_setup->GetLocalFieldManager(), true);

        //ic_em_field = new G4UniformElectricField(G4ThreeVector(0.0,1000.0*kilovolt/cm,0.0));
        //local_ic_field_mgr = new G4FieldManager(ic_em_field);
        //ion_cham_layer_log->SetFieldManager(local_ic_field_mgr, true);

        // Set as sensitive detector.
        ion_cham_layer_log->SetSensitiveDetector(ion_cham_layer_sd);

        // Make the ion chamber shell a sensitive volumes.
        if (! ion_cham_shell_sd)
          {
        	ion_cham_shell_sd = new SensitiveDet("/detector/ion_cham_shell");
            SDman->AddNewDetector(ion_cham_shell_sd);
          }
        ion_cham_shell_log->SetSensitiveDetector(ion_cham_shell_sd);

      }
    else
      {
        G4cout << "Excluding the ion chamber .........................................................." << G4endl;
      }


    if(Get_Sample_Flag() == "on")
      {

        G4cout << "Including the sample .................................................." << G4endl;

        // ------------------------------------- //
        // Sample foil.
        // ------------------------------------- //

        // Specify the dimensions of the sample foil.
        G4double sample_x_mm = 10.0 * mm;
		G4double sample_y_mm = 10.0 * mm;
		G4double sample_z_mm = 1.0 * mm;

		// Construct the sample foil.
        G4Box* sample = new G4Box("sample", sample_x_mm * 0.5, sample_y_mm * 0.5, sample_z_mm * 0.5);

		// Construct the logical volume for the sample foil.
		sample_log = new G4LogicalVolume(sample, sample_mat, "sample_log");

		// Construct the physical volume for the sample foil.
		sample_phys = new G4PVPlacement(0,              // Rotation.
										sample_pos,     // (x,y,z).
										sample_log,     // Logical volume.
										"sample_phys",  // Name.
										exp_hall_log,   // Mother volume.
										false,          // No boolean operations.
										0);             // Copy number.


		// Set the visualisation attributes of the ion chamber.
		G4VisAttributes* sample_visAtt = new G4VisAttributes(G4Colour(1.0,0.7,1.0));
        sample_log->SetVisAttributes(sample_visAtt);

        // Make the ion chamber a sensitive volume.
		if(! sample_sd)
			{
        	  sample_sd = new SensitiveDet("/detector/sample");
              SDman->AddNewDetector(sample_sd);
            }

        sample_log->SetSensitiveDetector(sample_sd);

      }
    else
    {
    	G4cout << "Excluding the sample .........................................................." << G4endl;
    }



  G4cout << "Detector construction complete" << G4endl;

  return exp_hall_phys;
}

void DetectorConstruction::Update_Geometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct_Geometry());
}

void DetectorConstruction::Set_Ion_Cham_Properties()
{
	// Map the strings to numbers as case switch works with ints.
	//map<const char*, int> ion_cham_map;
	map<G4String, unsigned int> ion_cham_map;
	ion_cham_map["LEFAC"] = 1;
	ion_cham_map["MEFAC"] = 2;

	// For now, hard code the ion chamber to always be the LEFAC.
	ion_cham_name = "LEFAC";

	// Get the integer value assigned to the ion chamber name.
	unsigned int input = ion_cham_map[ion_cham_name];

	// Pre-allocate a variable to which the name of the ion chamber fill gas will be assigned.
	G4String ion_cham_gas_mat_name;
	G4Material* ion_cham_gas_mat_ptr;
	G4double total_thick_mm;

	// Clear vectors.
	ion_cham_gas_mat_name_store.clear();
	ion_cham_gas_mat_store.clear();
	ion_cham_gas_thick_mm.clear();
	ion_cham_gas_shape_store.clear();
	ion_cham_gas_pos_store.clear();

	switch ( input ) {
	  case 1:
	    cout << "You selected LEFAC ...\n";

        // Specify the dimensions of the external section of the Al box that is the shell of the ion chamber.
	    ion_cham_ext_mm = G4ThreeVector(80.0*mm, 160.0*mm, 126.0*mm);

        // Specify the dimensions of the internal section of the Al box that is the shell of the ion chamber.
        // The shell thickness is 01.0*mm all the way around.
        ion_cham_wall_thick_mm = 10.0*mm;
        ion_cham_int_mm = G4ThreeVector(ion_cham_ext_mm.getX() - (2.0 * ion_cham_wall_thick_mm),
        		                        ion_cham_ext_mm.getY() - (2.0 * ion_cham_wall_thick_mm),
        		                        ion_cham_ext_mm.getZ() - (2.0 * ion_cham_wall_thick_mm));


        // Specify the dimensions of the sensitive volume of the ion chamber.
        ion_cham_sens_mm = G4ThreeVector(ion_cham_int_mm.getX() - correc_fac,
          		                         ion_cham_int_mm.getY() - correc_fac,
           		                         ion_cham_int_mm.getZ() - correc_fac);

	    // Store the gas material.
	    ion_cham_gas_mat_name = "Air";
        ion_cham_gas_mat_name_store.push_back(ion_cham_gas_mat_name);

        // Get the pointer to store.
	    ion_cham_gas_mat_ptr = G4Material::GetMaterial(ion_cham_gas_mat_name);

	    // Append the same material to all three layers.
	    ion_cham_gas_mat_store.push_back(ion_cham_gas_mat_ptr);
	    ion_cham_gas_mat_store.push_back(ion_cham_gas_mat_ptr);
	    ion_cham_gas_mat_store.push_back(ion_cham_gas_mat_ptr);

	    // Store the number of layers.
	    num_gas_layers = 3;
	    Set_Num_Gas_Layers(num_gas_layers);

        // Store the thicknesses of the layers.
        ion_cham_gas_thick_mm.push_back(43.0*mm);
        ion_cham_gas_thick_mm.push_back(20.0*mm);
        ion_cham_gas_thick_mm.push_back(43.0*mm);

        // Check that the number of layers is correct.
        assert (ion_cham_gas_thick_mm.size() == (unsigned int) num_gas_layers);

        // Check that the sum of the thicknesses fits in the volume allocated.
        total_thick_mm = std::accumulate(ion_cham_gas_thick_mm.begin(), ion_cham_gas_thick_mm.end(), 0);
        assert (total_thick_mm <= ion_cham_ext_mm.getZ());

        // Store the layer shapes.

        // Layer 1
        ion_cham_gas_shape_store.push_back(G4ThreeVector(ion_cham_sens_mm.getX() - correc_fac,
        		                                         ion_cham_sens_mm.getY() - correc_fac,
        		                                         ion_cham_gas_thick_mm[0] - correc_fac));

        // Layer 2
		ion_cham_gas_shape_store.push_back(G4ThreeVector(ion_cham_sens_mm.getX() - correc_fac,
				                                         ion_cham_sens_mm.getY() - correc_fac,
				                                         ion_cham_gas_thick_mm[1] - correc_fac ));
		// Layer 3
		ion_cham_gas_shape_store.push_back(G4ThreeVector(ion_cham_sens_mm.getX() - correc_fac,
				                                         ion_cham_sens_mm.getY() - correc_fac,
				                                         ion_cham_gas_thick_mm[2] - correc_fac));

        // Store the positions of the three layers.
		// These positions are absolute positions in 3D space and so must include the "ion_cham_pos" parameter.
        ion_cham_gas_pos_store.push_back(ion_cham_pos + G4ThreeVector(0.0, 0.0, -(0.5 * ion_cham_sens_mm.getZ()) + (0.5 * ion_cham_gas_shape_store[0].getZ()) + (0.5 * correc_fac)));
        ion_cham_gas_pos_store.push_back(ion_cham_pos + G4ThreeVector(0.0, 0.0, -(0.5 * ion_cham_sens_mm.getZ()) + (ion_cham_gas_shape_store[0].getZ()) + (0.5 * ion_cham_gas_shape_store[1].getZ()) + (1.0 * correc_fac)));
        ion_cham_gas_pos_store.push_back(ion_cham_pos + G4ThreeVector(0.0, 0.0, -(0.5 * ion_cham_sens_mm.getZ()) + (ion_cham_gas_shape_store[0].getZ()) + (ion_cham_gas_shape_store[1].getZ()) + (0.5 * ion_cham_gas_shape_store[2].getZ()) + (1.5 * correc_fac)));

        // Define the dimensions of the entrance window.
        ion_cham_ent_win_rad_mm = 6.25*mm;  // Radius of entrance window in mm.

        // Define the dimensions of the entrance window.
        ion_cham_exit_win_rad_mm = 20.0*mm; // Radius of exit window in mm.


	    break;
	  case 2:
	    cout << "You selected MEFAC ...\n";
	    break;
	  default:
	    cout<<"Error, bad input, quitting\n";
	    break;
	}

}

