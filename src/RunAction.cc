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
    // Store the base file name.
	file_base = val;

	// Store the file extension.
	file_ext = ".det";

	// Combine the base path / name and extension.
	data_file_name = file_base;
	data_file_name.insert(data_file_name.length(), file_ext);

	// Get new output file stream.
	data_out = new ofstream(data_file_name);

	// Store the file pointer.
	Store_File_Ptr(*data_out);

	// Store the file name.
	Store_Data_File_Name(data_file_name);
}

void RunAction::Set_Description(G4String val)
{
   description = val;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
	file_base = "/home/dimmockm/test";
	Set_File(file_base);
}

//----------------------------------------------------------------------------

void RunAction::EndOfRunAction(const G4Run*)
{
	if (file_base != ""){

		// Close the data file.
		data_out->close();

		// Delete the data file.
		delete data_out;
	}
}

