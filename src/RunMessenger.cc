#include "RunMessenger.hh"
#include "RunAction.hh"
#include "G4ios.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

RunMessenger::RunMessenger(RunAction * runinf)
:target(runinf)
{

  file_cmd = new G4UIcmdWithAString("/runinf/name",this);
  file_cmd->SetGuidance("Specify output file name");
  file_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

RunMessenger::~RunMessenger()
{
  delete file_cmd;
}

void RunMessenger::SetNewValue(G4UIcommand * command, G4String newValue)
{

  if( command == file_cmd )
    {
      target->Set_File(newValue);
    }
  G4cout << "Files saved as " << newValue << G4endl;

}

