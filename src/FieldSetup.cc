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
/// \file field/field03/src/FieldSetup.cc
/// \brief Implementation of the FieldSetup class
//
//
// $Id$
//
//  
//   Field Setup class implementation.
//

#include "FieldSetup.hh"
#include "FieldMessenger.hh"

#include "G4UniformMagField.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
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
#include "G4SystemOfUnits.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Constructors:

FieldSetup::FieldSetup()
  :  fChordFinder(0), fLocalChordFinder(0), fStepper(0)
{

  fEMField = new G4UniformElectricField(G4ThreeVector(0.0,100000.0*kilovolt/cm,0.0));

  fLocalEMField = new G4UniformElectricField(G4ThreeVector(0.0,100000.0*kilovolt/cm,0.0));

  fFieldMessenger = new FieldMessenger(this) ;
 
  fEquation = new G4EqMagElectricField(fEMField);
  fLocalEquation = new G4EqMagElectricField(fLocalEMField);
 
  fMinStep     = 0.25*mm ; // minimal step of 1 mm is default
  fStepperType = 4 ;      // ClassicalRK4 is default stepper

  fFieldManager = GetGlobalFieldManager();
  fLocalFieldManager = new G4FieldManager();

  UpdateField();
}

/////////////////////////////////////////////////////////////////////////////////

FieldSetup::FieldSetup(G4ThreeVector fieldVector)
{    
  fEMField = new G4UniformElectricField(fieldVector);
  fEquation = new G4EqMagElectricField(fEMField);
  fFieldManager = GetGlobalFieldManager();
  fFieldMessenger = new FieldMessenger(this) ;
  UpdateField();
}

////////////////////////////////////////////////////////////////////////////////

FieldSetup::~FieldSetup()
{
  if(fEMField) delete fEMField;
  if(fChordFinder)   delete fChordFinder;
  if(fStepper)       delete fStepper;
}

/////////////////////////////////////////////////////////////////////////////
//
// Update field
//

void FieldSetup::UpdateField()
{
  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

  fFieldManager->SetDetectorField(fEMField );
  fLocalFieldManager->SetDetectorField(fLocalEMField );

  if(fChordFinder) delete fChordFinder;
  if(fLocalChordFinder) delete fLocalChordFinder;

  fIntgrDriver = new G4MagInt_Driver(fMinStep,
                                       fStepper,
                                       fStepper->GetNumberOfVariables());

  fChordFinder = new G4ChordFinder(fIntgrDriver);

  fLocalIntgrDriver = new G4MagInt_Driver(fMinStep,
                                       fStepper,
                                       fStepper->GetNumberOfVariables());

  fLocalChordFinder = new G4ChordFinder(fLocalIntgrDriver);

  fFieldManager->SetChordFinder( fChordFinder );
  fLocalFieldManager->SetChordFinder( fLocalChordFinder );
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//

void FieldSetup::SetStepper()
{
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
}

/////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field to fieldValue along Z
//

void FieldSetup::SetFieldValue(G4double fieldStrength)
{
  G4ThreeVector fieldSetVec(0.0, 0.0, fieldStrength);
  this->SetFieldValue( fieldSetVec ); 
  //    *************

}

///////////////////////////////////////////////////////////////////////////////
//
// Set the value of the Global Field
//

void FieldSetup::SetFieldValue(G4ThreeVector fieldVector)
{
  if(fEMField) delete fEMField;

  if(fieldVector != G4ThreeVector(0.,0.,0.))
  { 
    fEMField = new  G4UniformElectricField(fieldVector);
  }
  else 
  {
    // If the new field's value is Zero, then 
    //  setting the pointer to zero ensures 
    //  that it is not used for propagation.
    fEMField = 0;
  }

  // Either  
  //   - UpdateField() to reset all (ChordFinder, Equation);
     // UpdateField();
  //   or simply update the field manager & equation of motion 
  //      with pointer to new field
  GetGlobalFieldManager()->SetDetectorField(fEMField);
  fEquation->SetFieldObj( fEMField );

}

////////////////////////////////////////////////////////////////////////////////
//
//  Utility method

G4FieldManager*  FieldSetup::GetGlobalFieldManager()
{
  return G4TransportationManager::GetTransportationManager()
                                ->GetFieldManager();
}


// In place of G4UniformField::GetConstantFieldValue ...
// 
G4ThreeVector FieldSetup::GetConstantFieldValue()
{
  static G4double fieldValue[6],  position[4]; 
  position[0] = position[1] = position[2] = position[3] = 0.0; 

  fEMField->GetFieldValue( position, fieldValue);
  G4ThreeVector fieldVec(fieldValue[0], fieldValue[1], fieldValue[2]); 

  return fieldVec;
}
