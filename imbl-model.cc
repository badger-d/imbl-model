//#############################################################################
// file name: silicon_wafers.cc
// This is the main file for the Si/HPD model. If you want to run this model
// on some older versions of Geant4, remove all of the "CLHEP::" identifiers 
// and add "using namespace CLHEP" at the top of the file.
// You might also need to swap any "G4VisExecutive" reference with "G4VisManager"
//
//
// Author: Toby Beveridge
//#############################################################################

#include <boost/cstdint.hpp>
#include <boost/date_time.hpp>
#include <boost/thread.hpp>

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh" 
#include <time.h>		
#include <sys/time.h>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif


//----------------------------------------------------------------------------

int main(int argc,char** argv) {
    const boost::system_time timeStarted = boost::get_system_time();

   //choose the Random engine
   CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

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

   //Initialise relevant classes ------------------------------------------
   // Run manager
   G4RunManager * runManager = new G4RunManager;

   //UserInitialization classes (mandatory)
   DetectorConstruction* detector = new DetectorConstruction;
   runManager->SetUserInitialization(detector);
   //PhysicsList* physics = new PhysicsList;
   //runManager->SetUserInitialization(physics);
   runManager->SetUserInitialization(new PhysicsList());
  
#ifdef G4VIS_USE
   //Visualization, if you choose to have it!
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();
#endif

  //UserAction classes
   PrimaryGeneratorAction* generator = new PrimaryGeneratorAction(detector);
   RunAction* run = new RunAction(detector,generator);

   runManager->SetUserAction(generator);
   runManager->SetUserAction(new EventAction(generator, run, detector));
   runManager->SetUserAction(run);

   //Initialize G4 kernel
//   runManager->Initialize();

   //get the pointer to the User Interface manager 
   G4UImanager * UI = G4UImanager::GetUIpointer();  

   //start a session (interactive or batch)------------------------------------------
   if(argc==1)
   // Define (G)UI terminal for interactive mode  
   { 
    	// G4UIterminal is a (dumb) terminal.
	G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
   session = new G4UIterminal(new G4UItcsh);      
#else
   session = new G4UIterminal();
#endif    
   //John change 26MAY09 --> don't default to using a vis file as it is bad in big runs
   //UI->ApplyCommand("/control/execute vis.mac");    
   session->SessionStart();
    	delete session;
   }
   else
   // Batch mode
   { 
	G4String command = "/control/execute ";
	G4String fileName = argv[1];
	UI->ApplyCommand(command+fileName);
   }
   //----------------------------------------------------------------------------------
   //tidy up:
#ifdef G4VIS_USE
   delete visManager;
#endif
   delete runManager;
   
   const boost::system_time timeEnded = boost::get_system_time();
   const uint64_t timeTakenMs = (timeEnded - timeStarted).total_milliseconds();

   G4cout << "Total time: " << timeTakenMs << " ms" << std::endl;

   return 0;
}

//----------------------------------------------------------------------------
