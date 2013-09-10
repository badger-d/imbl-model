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

  virtual void Set_File(G4String);

  virtual void Store_File_Ptr(std::ofstream &val){dataOut_ptr = &(val);};
  virtual G4String Get_File_Name(){return fileName;};
  virtual ofstream* Get_File_Ptr(){return dataOut_ptr;};

  private:
    DetectorConstruction* detector;
    PrimaryGeneratorAction* primary;
    RunMessenger* messenger;


    G4String fileBase;
    G4String fileExt;

    G4String fileName;
    ofstream* dataOut;
    ofstream* dataOut_ptr;


};

#endif





