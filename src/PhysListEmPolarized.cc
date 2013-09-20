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
// $Id: PhysListEmPolarized.cc,v 1.24 2009-11-15 22:10:03 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmPolarized.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4KleinNishinaModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

// Livermore replacements
#include "G4LivermorePolarizedPhotoElectricModel.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh"

#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LowEPComptonModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPolarized::PhysListEmPolarized(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPolarized::~PhysListEmPolarized()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmPolarized::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // Add Polarized EM Processes
  //
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      G4RayleighScattering* rs = new G4RayleighScattering;
      //rs->SetEmModel(new G4LivermoreRayleighModel());
      rs->SetEmModel(new G4LivermorePolarizedRayleighModel());
      ph->RegisterProcess(rs, particle);

      G4PhotoElectricEffect* pe = new G4PhotoElectricEffect;
      //pe->SetEmModel(new G4LivermorePhotoElectricModel());
      pe->SetEmModel(new G4LivermorePolarizedPhotoElectricModel());
      ph->RegisterProcess(pe, particle);            

      G4ComptonScattering* cs   = new G4ComptonScattering;
      //cs->SetEmModel(new G4LowEPComptonModel());
      cs->SetEmModel(new G4LivermorePolarizedComptonModel());
      ph->RegisterProcess(cs, particle);

      G4GammaConversion* conv = new G4GammaConversion();
      //conv->SetEmModel(new G4LivermoreGammaConversionModel());
      conv->SetEmModel(new G4LivermorePolarizedGammaConversionModel());
      ph->RegisterProcess(conv, particle);
     
    } else if (particleName == "e-") {

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);      
      ph->RegisterProcess(eIoni,                   particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);      
	    
    } else if (particleName == "e+") {

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);      
      ph->RegisterProcess(eIoni,                     particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),   particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);
                  
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.1, 50*um);      
      ph->RegisterProcess(muIoni,                   particle);
      ph->RegisterProcess(new G4MuBremsstrahlung(), particle);
      ph->RegisterProcess(new G4MuPairProduction(), particle);
                   
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {
  
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.1, 20*um);
      ph->RegisterProcess(hIoni,                   particle);
      ph->RegisterProcess(new G4hBremsstrahlung(), particle);
      ph->RegisterProcess(new G4hPairProduction(), particle);            
     
    } else if( particleName == "alpha" || 
	       particleName == "He3"    ) {

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(ionIoni,                 particle);
      ph->RegisterProcess(new G4NuclearStopping(), particle);      
            
    } else if( particleName == "GenericIon" ) {
 
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(ionIoni,                 particle);
      ph->RegisterProcess(new G4NuclearStopping(), particle);                   
      
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
	       
      //all others charged particles except geantino
      ph->RegisterProcess(new G4hIonisation(),    particle);
    }
  }
    
  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  de->SetFluo(true);
  de->SetAuger(true);   
  de->SetPIXE(true);  
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

