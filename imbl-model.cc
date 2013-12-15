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

#include <boost/cstdint.hpp>
#include <boost/date_time.hpp>
#include <boost/thread.hpp>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv) 
{

	const boost::system_time timeStarted = boost::get_system_time();

	// Choose the random engine.
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

	// Generate seed from system file.
	uint64_t seed = 0;
	std::ifstream randomFile;

	randomFile.exceptions(~std::ios::goodbit);
	randomFile.open("/dev/random", std::ios::binary);
	randomFile.read((char *) & seed, sizeof(seed));
	randomFile.close();

	G4cout << "Random seed is " << seed << G4endl;

	if (argc == 3) {
		// second argument is assumed to be random seed
		seed = atoi(argv[2]);
		G4cout << "Input random seed is " << seed << G4endl;
	}

	CLHEP::HepRandom::setTheSeed(seed);
	CLHEP::HepRandom::showEngineStatus();

	//CLHEP::HepRandom::saveEngineStatus("./currentEvent-MRD.rndm");
	//CLHEP::HepRandom::restoreEngineStatus("./currentEvent-MRD.rndm");

	// Initialise relevant classes ------------------------------------------

	// Init. run manager.
	G4RunManager* runManager = new G4RunManager;

	// Init. detector construction.
	DetectorConstruction* detector = new DetectorConstruction;

	// Pass the detector into the run.
	runManager->SetUserInitialization(detector);

	// Init. physics list.
	PhysicsList* physics = new PhysicsList();

	// Load physics list.
	runManager->SetUserInitialization(physics);

	// Init. primary generator.
	PrimaryGeneratorAction* generator = new PrimaryGeneratorAction(detector);

	// Load primary generator.
	runManager->SetUserAction(generator);

	// Init. run action.
	RunAction* runAction = new RunAction(detector, generator);

	// Load run action.
	runManager->SetUserAction(runAction);

	// Init. event action.
	EventAction* eventAction = new EventAction(runAction, detector);

	// Load event action.
	runManager->SetUserAction(eventAction);

	// G4 kernel, physics tables init. in macro so I've commented out the line.
	//runManager->Initialize();

#ifdef G4VIS_USE

	// visualization manager

	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();

#endif 

	// Get the pointer to the User Interface manager

	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	if (argc!=1)   // batch mode
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);
	}
	else
	{  // interactive mode : define UI session
#ifdef G4UI_USE
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
		ui->SessionStart();
		delete ui;
#endif
	}

	// job termination

#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;

	const boost::system_time timeEnded = boost::get_system_time();
	const uint64_t timeTakenMs = (timeEnded - timeStarted).total_milliseconds();

	G4cout << "Total time: " << timeTakenMs << " ms" << std::endl;


	return 0;
}
