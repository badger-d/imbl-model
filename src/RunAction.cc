#include "RunAction.hh"
#include "RunMessenger.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <time.h>
#include <boost/date_time.hpp>
#include <boost/thread.hpp>

RunAction::RunAction(DetectorConstruction* DC, PrimaryGeneratorAction* PG)
{
   messenger = new RunMessenger(this);
   detector = DC;
   primary = PG;
}

RunAction::~RunAction()
{
	// Delete the run messenger.
	delete messenger;
}

void RunAction::Set_File(G4String val)
{
    fileBase = val;

	// Prep the data file.
	fileExt = ".det";
	fileName = fileBase;
	fileName.insert(fileName.length(), fileExt);
	dataOut = new ofstream(fileName);
	Store_File_Ptr(*dataOut);
}

void RunAction::BeginOfRunAction(const G4Run*)
{

}

//----------------------------------------------------------------------------

void RunAction::EndOfRunAction(const G4Run*)
{
	if (fileBase != ""){

		// Close the data file.
		dataOut->close();

		// Delete the data file.
		delete dataOut;
	}
}

