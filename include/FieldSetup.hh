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
/// \file field/field03/include/FieldSetup.hh
/// \brief Definition of the FieldSetup class
//
// $Id$
//

#ifndef FieldSetup_H
#define FieldSetup_H

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"

class FieldMessenger;
class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4MagInt_Driver;

///  A class for setting up the Magnetic Field
///
///  It also creates the necessary classes to control accuracy of propagation.
///  In this example
///    - There is a global field for most of the setup;
///    - A local  field overides it for some volume(s) and it assumed to be 
///      uniform.

class FieldSetup
{

public:

  FieldSetup();               //  A zero field
  FieldSetup(G4ThreeVector);  //  The value of the field
  ~FieldSetup();
      
  void SetStepperType(G4int i) { fStepperType = i; }

  void SetStepper();

  void SetMinStep(G4double s) { fMinStep = s; }

  void UpdateField();

  void SetFieldValue(G4ThreeVector fieldVector);
  void SetFieldValue(G4double      fieldValue);
  G4ThreeVector GetConstantFieldValue();
  G4FieldManager*  GetLocalFieldManager() { return fLocalFieldManager;}

protected:

  G4FieldManager*         GetGlobalFieldManager();
    // Returns the global Field Manager

  G4FieldManager*         fFieldManager;
  G4FieldManager*         fLocalFieldManager;

  G4ChordFinder*          fChordFinder;
  G4ChordFinder*          fLocalChordFinder;

  G4EqMagElectricField*   fEquation;
  G4EqMagElectricField*   fLocalEquation;

  G4UniformElectricField* fEMField;
  G4UniformElectricField* fLocalEMField;

  G4MagInt_Driver*        fIntgrDriver;
  G4MagInt_Driver*        fLocalIntgrDriver;

  G4MagIntegratorStepper* fStepper;
  G4MagIntegratorStepper* fLocalStepper;
  G4int                   fStepperType;

  G4double                fMinStep;
 
  FieldMessenger*      fFieldMessenger;

};

#endif
