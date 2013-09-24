#ifndef RunMessenger_h
#define RunMessenger_h 1

class RunAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class RunMessenger: public G4UImessenger
{
  public:
    RunMessenger(RunAction* runinf);
    ~RunMessenger();

    virtual void SetNewValue(G4UIcommand * command,G4String newValues);

  private:
    RunAction* target;
    //commands
    G4UIcmdWithAString*     file1Cmd;
    G4UIcmdWithAString*     descriptionCmd;
};

#endif
