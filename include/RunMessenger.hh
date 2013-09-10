// RunMessenger.hh,v 1.3 2002/12/13 11:34:28 gunter Exp $
// --------------------------------------------------------------
//
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

    virtual void SetNewValue(G4UIcommand * command, G4String newValues);

  private:
    RunAction* target;

    G4UIcmdWithAString*     file_cmd;
};

#endif

