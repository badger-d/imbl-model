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
// User Field class implementation.
//

#include "ElectricFieldSetup.hh"
#include "FieldMessenger.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//  Constructors:

ElectricFieldSetup::ElectricFieldSetup()
: fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fEMfield(0),
   fElFieldValue(),
   fStepper(0),
   fIntgrDriver(0),
   fStepperType(4),    // ClassicalRK4 -- the default stepper
   fMinStep(0.010*mm)  // minimal step of 10 microns
 {

	  fEMfield = new G4UniformElectricField(
	                   G4ThreeVector(0.0,100000.0*kilovolt/cm,0.0));
	  fEquation = new G4EqMagElectricField(fEMfield);

	  fFieldManager = GetGlobalFieldManager();
	  fFieldMessenger = new FieldMessenger(this) ;
	  UpdateField();


//  fChordFinder = 0;
//  //fLocalChordFinder = 0;
//  fStepper = 0;
//
//  fEMfield = new G4UniformElectricField(0.0,1.0*kilovolt/cm,0.0);
//  //fLocalEMfield = new G4UniformElectricField(0.0,1.0*kilovolt/cm,0.0);
//
//  fFieldMessenger = new FieldMessenger(this);
//
//  fEquation = new G4EqMagElectricField(fEMfield);
//  //fLocalEquation = new G4EqMagElectricField(fLocalEMfield);
//
//  fMinStep     = 0.25*mm ; // minimal step of 1 mm is default
//  fStepperType = 4 ;      // ClassicalRK4 is default stepper
//
//  fFieldManager = GetGlobalFieldManager();
//  //fLocalFieldManager = new G4FieldManager();
//
//  UpdateField();










}

ElectricFieldSetup::ElectricFieldSetup(G4ThreeVector fieldVector)
: fFieldManager(0),
   fChordFinder(0),
   fEquation(0),
   fEMfield(0),
   fElFieldValue(),
   fStepper(0),
   fIntgrDriver(0),
   fStepperType(4),    // ClassicalRK4 -- the default stepper
   fMinStep(0.010*mm)  // minimal step of 10 microns
{    
	  fEMfield = new G4UniformElectricField(fieldVector);
	  // GetGlobalFieldManager()->CreateChordFinder(this);
	  fEquation = new G4EqMagElectricField(fEMfield);

	  fFieldManager = GetGlobalFieldManager();
	  fFieldMessenger = new FieldMessenger(this) ;
	  UpdateField();

//	fEMfield = new G4UniformElectricField(fieldVector);
//	///////// GetGlobalFieldManager()->CreateChordFinder(this);
//	fEquation = new G4EqMagElectricField(fEMfield);
//	fFieldMessenger = new FieldMessenger(this);
//	fFieldManager = GetGlobalFieldManager();
//	UpdateField();








}

ElectricFieldSetup::~ElectricFieldSetup()
{
if (fChordFinder) delete fChordFinder;
if (fStepper)     delete fStepper;
if (fEquation)    delete fEquation;
if (fEMfield)     delete fEMfield;

//	  if(fEMfield) delete fEMfield;
//	  if(fChordFinder)   delete fChordFinder;
//	  if(fIntgrDriver) delete fIntgrDriver;
//	  if(fStepper)       delete fStepper;








}

void ElectricFieldSetup::UpdateField()
{




	// Register this field to 'global' Field Manager and
	// Create Stepper and Chord Finder with predefined type, minstep (resp.)

	  SetStepper();

	  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	  fFieldManager->SetDetectorField(fEMfield );

	  if (fChordFinder) delete fChordFinder;
	  // fChordFinder = new G4ChordFinder( fEMfield, fMinStep, fStepper);

	  fIntgrDriver = new G4MagInt_Driver(fMinStep,
	                                     fStepper,
	                                     fStepper->GetNumberOfVariables());

	  fChordFinder = new G4ChordFinder(fIntgrDriver);

	  fFieldManager->SetChordFinder(fChordFinder);











	//// Register this field to 'global' Field Manager and
//// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//
//	  SetStepper();
//	  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
//
//	  fFieldManager->SetDetectorField(fEMfield );
//	  //fLocalFieldManager->SetDetectorField(fLocalEMfield );
//
//	  if(fChordFinder) delete fChordFinder;
//	 // if(fLocalChordFinder) delete fLocalChordFinder;
//      G4cout << " MOFO " << fStepper << " " << fStepper->GetNumberOfVariables() << " " << fStepperType << G4endl;
//	  fIntgrDriver = new G4MagInt_Driver(fMinStep,
//	                                     fStepper,
//	                                     fStepper->GetNumberOfVariables());
//
//      //fLocalIntgrDriver = new G4MagInt_Driver(fMinStep,
//	  //                                        fLocalStepper,
//	  //                                        fLocalStepper->GetNumberOfVariables());
//
//      fChordFinder = new G4ChordFinder(fIntgrDriver);
//
//	  //fLocalChordFinder = new G4ChordFinder(fLocalIntgrDriver);
//
//	  fFieldManager->SetChordFinder( fChordFinder );
//	  //fLocalFieldManager->SetChordFinder( fLocalChordFinder );









}

void ElectricFieldSetup::SetStepper()
{






	// Set stepper according to the stepper type

	  G4int nvar = 8;

	  if (fStepper) delete fStepper;

	  switch ( fStepperType )
	  {
	    case 0:
	      fStepper = new G4ExplicitEuler( fEquation, nvar );
	      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;
	      break;
	    case 1:
	      fStepper = new G4ImplicitEuler( fEquation, nvar );
	      G4cout<<"G4ImplicitEuler is called"<<G4endl;
	      break;
	    case 2:
	      fStepper = new G4SimpleRunge( fEquation, nvar );
	      G4cout<<"G4SimpleRunge is called"<<G4endl;
	      break;
	    case 3:
	      fStepper = new G4SimpleHeum( fEquation, nvar );
	      G4cout<<"G4SimpleHeum is called"<<G4endl;
	      break;
	    case 4:
	      fStepper = new G4ClassicalRK4( fEquation, nvar );
	      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
	      break;
	    case 5:
	      fStepper = new G4CashKarpRKF45( fEquation, nvar );
	      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
	      break;
	    case 6:
	      fStepper = 0; // new G4RKG3_Stepper( fEquation, nvar );
	      G4cout<<"G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
	      break;
	    case 7:
	      fStepper = 0; // new G4HelixExplicitEuler( fEquation );
	      G4cout<<"G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
	      break;
	    case 8:
	      fStepper = 0; // new G4HelixImplicitEuler( fEquation );
	      G4cout<<"G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
	      break;
	    case 9:
	      fStepper = 0; // new G4HelixSimpleRunge( fEquation );
	      G4cout<<"G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
	      break;
	    default: fStepper = 0;

	  }









//	  G4int nvar = 8;
//
//	  if(fStepper) delete fStepper;
//      switch ( fStepperType )
//	  {
//	  case 0:
//	      fStepper = new G4ExplicitEuler( fEquation, nvar );
//	      //fLocalStepper = new G4ExplicitEuler( fLocalEquation, nvar );
//	      G4cout<<"G4ExplicitEuler is calledS"<<G4endl;
//	      break;
//	    case 1:
//	      fStepper = new G4ImplicitEuler( fEquation, nvar );
//	      //fLocalStepper = new G4ImplicitEuler( fLocalEquation, nvar );
//	      G4cout<<"G4ImplicitEuler is called"<<G4endl;
//	      break;
//	    case 2:
//	      fStepper = new G4SimpleRunge( fEquation, nvar );
//	      //fLocalStepper = new G4SimpleRunge( fLocalEquation, nvar );
//	      G4cout<<"G4SimpleRunge is called"<<G4endl;
//	      break;
//	    case 3:
//	      fStepper = new G4SimpleHeum( fEquation, nvar );
//	      //fLocalStepper = new G4SimpleHeum( fLocalEquation, nvar );
//	      G4cout<<"G4SimpleHeum is called"<<G4endl;
//	      break;
//	    case 4:
//	      fStepper = new G4ClassicalRK4( fEquation, nvar );
//	      //fLocalStepper = new G4ClassicalRK4( fLocalEquation, nvar );
//	      G4cout<<"G4ClassicalRK4 (default) is called"<<G4endl;
//	      break;
//	    case 5:
//	      fStepper = new G4CashKarpRKF45( fEquation, nvar );
//	      //fLocalStepper = new G4CashKarpRKF45( fLocalEquation, nvar );
//	      G4cout<<"G4CashKarpRKF45 is called"<<G4endl;
//	      break;
//	    case 6:
//	      fStepper = 0; // new G4RKG3_Stepper( fEquation, nvar );
//	      //fLocalStepper = 0;
//	      G4cout<<"G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
//	      break;
//	    case 7:
//	      fStepper = 0; // new G4HelixExplicitEuler( fEquation );
//	      //fLocalStepper = 0;
//	      G4cout<<"G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
//	      break;
//	    case 8:
//	      fStepper = 0; // new G4HelixImplicitEuler( fEquation );
//	      //fLocalStepper = 0;
//	      G4cout<<"G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
//	      break;
//	    case 9:
//	      fStepper = 0; // new G4HelixSimpleRunge( fEquation );
//	      //fLocalStepper = 0;
//	      G4cout<<"G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
//	      break;
//
//		default: fStepper = 0;//, fLocalStepper = 0;
//	  }




}



void ElectricFieldSetup::SetFieldValue(G4double fieldStrength)
{




	// Set the value of the Global Field to fieldValue along Z

	  G4ThreeVector fieldVector( 0.0, 0.0, fieldStrength );

	  SetFieldValue( fieldVector );








	//
//	  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
//	  this->SetFieldValue( fieldSetVec );
//










}

void ElectricFieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{





	// Set the value of the Global Field value to fieldVector

	  // Find the Field Manager for the global field
	  G4FieldManager* fieldMgr= GetGlobalFieldManager();

	  if (fieldVector != G4ThreeVector(0.,0.,0.))
	  {
	    if (fEMfield) delete fEMfield;
	    fEMfield = new G4UniformElectricField(fieldVector);

	    fEquation->SetFieldObj(fEMfield);  // must now point to the new field

	    // UpdateField();

	    fieldMgr->SetDetectorField(fEMfield);
	  }
	  else
	  {
	    // If the new field's value is Zero, then it is best to
	    //  insure that it is not used for propagation.
	    if (fEMfield) delete fEMfield;
	    fEMfield = 0;
	    fEquation->SetFieldObj(fEMfield);   // As a double check ...
	    fieldMgr->SetDetectorField(fEMfield);
	  }















//
//
//	  if(fEMfield) delete fEMfield;
//
//	  if(fieldVector != G4ThreeVector(0.,0.,0.))
//	  {
//	    fEMfield = new  G4UniformElectricField(fieldVector);
//	  }
//	  else
//	  {
//	    // If the new field's value is Zero, then
//	    //  setting the pointer to zero ensures
//	    //  that it is not used for propagation.
//	    fEMfield = 0;
//	  }
//
//	  // Either
//	  //   - UpdateField() to reset all (ChordFinder, Equation);
//	     // UpdateField();
//	  //   or simply update the field manager & equation of motion
//	  //      with pointer to new field
//	  GetGlobalFieldManager()->SetDetectorField(fEMfield);
//	  fEquation->SetFieldObj( fEMfield );
//
//










}

G4FieldManager*  ElectricFieldSetup::GetGlobalFieldManager()
{



	  return G4TransportationManager::GetTransportationManager()
	           ->GetFieldManager();





//  return G4TransportationManager::GetTransportationManager()
//                          ->GetFieldManager();






}

//G4ThreeVector ElectricFieldSetup::GetConstantFieldValue()
//{
//  static G4double fieldValue[6],  position[4];
//  position[0] = position[1] = position[2] = position[3] = 0.0;
//
//  fEMfield->GetFieldValue( position, fieldValue);
//  G4ThreeVector fieldVec(fieldValue[0], fieldValue[1], fieldValue[2]);
//
//  return fieldVec;
//}

