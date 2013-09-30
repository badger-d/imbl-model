#include "RunMessenger.hh"
#include "RunAction.hh"
#include "G4ios.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

RunMessenger::RunMessenger(RunAction * runinf)
:target(runinf)
{

  file1Cmd = new G4UIcmdWithAString("/runinf/fname_conv",this);
  file1Cmd->SetGuidance("Specify output filename convention");
  file1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  descriptionCmd = new G4UIcmdWithAString("/runinf/description",this);
  descriptionCmd->SetGuidance("Specify description of this experiment");
  descriptionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

RunMessenger::~RunMessenger()
{

  delete file1Cmd;
  delete descriptionCmd;

}

void RunMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command == descriptionCmd )
    {
      target->Set_Description(newValue);
      G4cout << "Experiment description " << newValue << G4endl;
    }

  if( command == file1Cmd )
    {
      target->Set_File(newValue);
      G4cout << "Files saved as " << newValue << G4endl;
    }

}
