#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun)
:Action(Gun)
{

	// Set the particle energy
	energy_cmd = new G4UIcmdWithADoubleAndUnit("/gun/inpenergy", this);
	energy_cmd->SetGuidance("Select particle energy.");
	energy_cmd->SetParameterName("Size", false);
	energy_cmd->SetRange("Size>=0.");
	energy_cmd->SetDefaultValue(10.);
	energy_cmd->SetUnitCategory("Energy");
	energy_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	// Set the emission direction distribution.
	isotropy_cmd = new G4UIcmdWithAString("/gun/isotropy",this);
	isotropy_cmd->SetGuidance("Select the emission direction for barticles.");
	isotropy_cmd->SetGuidance("Choice : twopi(default), fourpi, forward");
	isotropy_cmd->SetParameterName("choice",true);
	isotropy_cmd->SetDefaultValue("twopi");
	isotropy_cmd->SetCandidates("twopi fourpi forward");
	isotropy_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	// Set the emission location.
	origin_cmd = new G4UIcmdWith3VectorAndUnit("/gun/origin", this);
	origin_cmd->SetGuidance("Select particle origin.");
	origin_cmd->SetParameterName("ox","oy","oz", true, true);
	origin_cmd->SetDefaultUnit("mm");

	// Set the diameter at the point of of the emission.
	diameter_cmd = new G4UIcmdWithADoubleAndUnit("/gun/diameter", this);
	diameter_cmd->SetGuidance("Select width of source that refers to diameter (disc) or side length (square).");
	diameter_cmd->SetParameterName("Size",false);
	diameter_cmd->SetRange("Size>=0.");
	diameter_cmd->SetUnitCategory("Length");
	diameter_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	// Set the emission position distribution.
	distribution_cmd = new G4UIcmdWithAString("/gun/distribution",this);
	distribution_cmd->SetGuidance("Shoot particles in random directions");
	distribution_cmd->SetGuidance("Choice : point, disc, square");
	distribution_cmd->SetParameterName("choice",true);
	distribution_cmd->SetDefaultValue("point");
	distribution_cmd->SetCandidates("point disc square");
	distribution_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete energy_cmd;
  delete isotropy_cmd;
  delete origin_cmd;
  delete diameter_cmd;
  delete distribution_cmd;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{

  if(command == energy_cmd)
    {
      Action->Set_Energy(energy_cmd->GetNewDoubleValue(newValue));
      G4cout << "Particle (beam) energy changed to " << newValue << " ............................" << G4endl;
    }

  if( command == isotropy_cmd )
      {
        Action->Set_Isotropy(newValue);
        G4cout << "Source emission vector isotropy changed to " << newValue << " ............................" << G4endl;
      }

  if(command == origin_cmd)
    {
      Action->Set_Origin(origin_cmd->GetNew3VectorValue(newValue));
      G4cout << "Emission position changed to " << newValue << " ............................" << G4endl;
    }

  if(command == diameter_cmd)
    {
      Action->Set_Diameter(diameter_cmd->GetNewDoubleValue(newValue));
      G4cout << "Beam diameter changed to " << newValue << " ............................" << G4endl;
    }

  if( command == distribution_cmd )
      {
        Action->Set_Distribution(newValue);
        G4cout << "Source emission location distribution changed to " << newValue << " ............................" << G4endl;
      }
}
