#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
:Action(Gun)
{

  energy_cmd = new G4UIcmdWithADoubleAndUnit("/gun/inpenergy", this);
  energy_cmd->SetGuidance("select particle energy");
  energy_cmd->SetParameterName("Size", false);
  energy_cmd->SetRange("Size>=0.");
  energy_cmd->SetDefaultValue(10.);
  energy_cmd->SetUnitCategory("Energy");
  energy_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  isotropy_cmd = new G4UIcmdWithAString("/gun/isotropy",this);
  isotropy_cmd->SetGuidance("Shoot particles in random directions");
  isotropy_cmd->SetGuidance("Choice : twopi(default), fourpi, fracpi");
  isotropy_cmd->SetParameterName("choice",true);
  isotropy_cmd->SetDefaultValue("twopi");
  isotropy_cmd->SetCandidates("twopi fourpi fracpi");
  isotropy_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  origin_cmd = new G4UIcmdWith3VectorAndUnit("/gun/origin", this);
  origin_cmd->SetGuidance("select particle origin");
  origin_cmd->SetParameterName("ox","oy","oz", true, true);
  origin_cmd->SetDefaultUnit("mm");

}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete energy_cmd;
  delete isotropy_cmd;
  delete origin_cmd;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{

  if(command == energy_cmd)
    {
      Action->Set_Energy(energy_cmd->GetNewDoubleValue(newValue));
      G4cout << "Particle energy changed to " << newValue << " ............................" << G4endl;
    }

  if( command == isotropy_cmd )
      {
        Action->Set_Isotropy(newValue);
        G4cout << "Source shape changed to " << newValue << " ............................" << G4endl;
      }

  if(command == origin_cmd)
    {
      Action->Set_Origin(origin_cmd->GetNew3VectorValue(newValue));
      G4cout << "Particle position changed to " << newValue << " ............................" << G4endl;
    }
}
