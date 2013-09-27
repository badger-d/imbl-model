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
  distribution = "point";
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

	// Specify the particle momentum.
	G4ThreeVector momentum = G4ThreeVector(1.0,0.0,0.0);  // Define momentum as a vector

	// Get the emission location from the desired distribution e.g. disc.
	// This will be centered around the origin.
	G4ThreeVector location = Distribution(distribution);   // D

	// Add the origin to the locaiton as an offset.
	location += origin;

	// Set the emission location to the correct value.
	particleGun->SetParticlePosition(location);

    // Specify variables for spherical coordinates of emission vector.
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

G4ThreeVector PrimaryGeneratorAction::Distribution(G4String shape)
{

	// Init. some vars.
	G4ThreeVector location = G4ThreeVector(0,0,0);
	G4double temp_r, temp_theta, temp_x, temp_y;

	// Set the location for a point source.
	if (shape == "point"){

		location.set(0.0, 0.0, 0.0);
		Set_Diameter(0.0);

	}	else if (shape == "disc") {

		// Sample the location for a disc source.

		temp_theta = 2 * G4UniformRand() * pi;
		temp_r = sqrt(G4UniformRand()) * diameter * 0.5;
		temp_x = (temp_r * (cos(temp_theta)));
		temp_y = (temp_r * (sin(temp_theta)));

		//G4cout << src_choice << G4endl;
		location.set(temp_x * mm, temp_y * mm, 0.0 * mm);

	} else if (shape == "square"){

		// Sample the location for a square source.

		temp_x = (G4UniformRand()-0.5) * diameter;
		temp_y = (G4UniformRand()-0.5) * diameter;

		location.set(temp_x, temp_y, 0.0);

	}	else {

		//Error
		exit(1);

	}

	//location.rotateZ(pi/2); // USE WITH PHANTOM
	return location;
}
