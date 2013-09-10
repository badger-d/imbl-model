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
  isotropy = "onepi";
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
	  theta = acos(2*(G4UniformRand())-1);
	  phi = 2*((pi)*(G4UniformRand()));
  }

  else if(isotropy == "twopi")
  {
	  theta = acos(2*(G4UniformRand())-1);
	  phi = pi + ((pi)*(G4UniformRand()));
  }

  else if(isotropy == "fracpi")
  {
	  // Generate vector in +/- 45deg squared cone from emission location
	  G4double const frac  = 0.1 * pi;
	  G4double const theta_mean = 0.5*pi;
	  G4double const phi_mean   = 1.5*pi;

	  theta = (G4UniformRand()-0.5) * frac + theta_mean;
	  phi   = (G4UniformRand()-0.5) * frac + phi_mean;
  }

  momentum.setRThetaPhi(r,theta,phi);
  momentum.rotateZ(pi/2.0);

  particleGun->SetParticleMomentumDirection(momentum);

  // Randomize polarization.
  G4ThreeVector polarization = momentum.orthogonal();
  (polarization.rotate(2*pi*G4UniformRand(), momentum)).unit();

  particleGun->SetParticlePolarization(polarization);
  particleGun->SetParticleMomentumDirection(momentum);
  particleGun->GeneratePrimaryVertex(anEvent);
}

