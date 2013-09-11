#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorMessenger.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

#include <sys/time.h>
#include <iostream>
#include <fstream>


PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC)
:detector(myDC)
{
  // Set initial values
  isotropy = "forward";
  energy = 1. * keV;
  origin = G4ThreeVector(0.,0.,0.);

  particleGun = new G4ParticleGun(1);
  messenger = new PrimaryGeneratorMessenger(this);

  // Default particle
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("gamma");

  // Gun specifications
  particleGun->SetParticlePosition(origin);
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleEnergy(energy);

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete messenger;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun->SetParticleEnergy(energy);
  particleGun->SetParticlePosition(origin);

  G4ThreeVector momentum = G4ThreeVector(1.0,0.0,0.0);  // Define momentum as a vector
  particleGun->SetParticlePosition(origin);

  G4double r = 1.0;
  G4double theta = 0.0;
  G4double phi = 0.0;


  if(isotropy == "fourpi")
  {
	  // Sample random emission direction in 4pi.
	  theta = acos(2*(G4UniformRand())-1);
	  phi = 2*((pi)*(G4UniformRand()));
  }

  else if(isotropy == "twopi")
  {
	  // Sample random emission direction in 2pi.
	  theta = acos(2*(G4UniformRand())-1);
	  phi = pi + ((pi)*(G4UniformRand()));
  }

  else if(isotropy == "forward")
  {
      // Emit particle in forward direction (-z to + z).
	  theta = acos(2*(0.5));
	  phi = pi;
  }

  momentum.setRThetaPhi(r,theta,phi);

  particleGun->SetParticleMomentumDirection(momentum);

  // Randomize polarization.
  G4ThreeVector polarization = momentum.orthogonal();
  (polarization.rotate(2*pi*G4UniformRand(), momentum)).unit();

  particleGun->SetParticlePolarization(polarization);
  particleGun->SetParticleMomentumDirection(momentum);
  particleGun->GeneratePrimaryVertex(anEvent);
}

