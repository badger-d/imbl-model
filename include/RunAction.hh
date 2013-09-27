#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
//----------------------------------------------------------------------------
using namespace std;
class G4Run;
class RunMessenger;
class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*,PrimaryGeneratorAction*);
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  // Store the pointer to the header file.
  virtual void Store_Head_File_Name(G4String val){stored_head_fname = val;};

  virtual G4String Get_Head_File_Name(){return stored_head_fname;};

  virtual G4String Get_Data_File_Name(){return stored_data_fname;};

  virtual  void Store_Experiment_Name(G4String val){stored_exp_name = val;};

  virtual  G4String Get_Experiment_Name(){return stored_exp_name;};

  virtual void Set_File(G4String);
  virtual void Set_Description(G4String);

  // Store the pointer to the output file.
  virtual void Store_File_Ptr(std::ofstream &val){data_out_ptr = &(val);};

  // Get the file name.
  virtual G4String Get_File_Name(){return file_name;};

  // Store the data-file name.
  virtual void Store_Data_File_Name(G4String val){stored_data_file_name = val;};

  // Get the output file pointer.
  virtual ofstream* Get_File_Ptr(){return data_out_ptr;}

  private:
    DetectorConstruction* detector;
    PrimaryGeneratorAction* primary;
    RunMessenger* messenger;

    G4String stored_head_fname;
    G4String stored_data_fname;
    G4String stored_exp_name;

    // Pointer to the output file pointer.
    ofstream* data_out_ptr;

    // Pointer to the output file.
    ofstream* data_out;

    G4String data_file_name;

    G4String stored_data_file_name;

    G4String file_base;
    G4String file_ext;

    G4String file_name;

    G4String description;


};

#endif





