#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4Material.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* myDet)
:myDetector(myDet)
{

  detector_cmd = new G4UIcmdWithAString("/geometry/detector",this);
  detector_cmd->SetGuidance("Selects inclusion of insulated box");
  detector_cmd->SetGuidance("Choice : on(default), off");
  detector_cmd->SetParameterName("choice",true);
  detector_cmd->SetDefaultValue("on");
  detector_cmd->SetCandidates("on off");
  detector_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sample_cmd = new G4UIcmdWithAString("/geometry/sample",this);
  sample_cmd->SetGuidance("Selects inclusion of sample");
  sample_cmd->SetGuidance("Choice : on(default), off");
  sample_cmd->SetParameterName("choice",true);
  sample_cmd->SetDefaultValue("on");
  sample_cmd->SetCandidates("on off");
  sample_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sampleX_cmd = new G4UIcmdWithADoubleAndUnit("/geometry/sample_x", this);
  sampleX_cmd->SetGuidance("Set depth of PD.");
  sampleX_cmd->SetParameterName("Size",false);
  sampleX_cmd->SetRange("Size>=0.");
  sampleX_cmd->SetUnitCategory("Length");
  sampleX_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  update_cmd = new G4UIcmdWithoutParameter("/geometry/update",this);
  update_cmd->SetGuidance("Update detector geometry.");
  update_cmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  update_cmd->SetGuidance("if you changed geometrical value(s).");
  update_cmd->AvailableForStates(G4State_Idle);

}


DetectorMessenger::~DetectorMessenger()
{

  delete detector_cmd;
  delete sample_cmd;
  delete sampleX_cmd;
  delete update_cmd;

}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{


  if(command == detector_cmd)
    {
      myDetector->Set_Ion_Cham_Flag(newValue);
      G4cout << "### Ion chamber included." << newValue << G4endl;
    }

  if(command == sample_cmd)
    {
	  myDetector->Set_Sample_Flag(newValue);
	  G4cout << "### Sample included." << newValue << G4endl;
    }

  if(command == sampleX_cmd)
    {
      myDetector->Set_Sample_X(sampleX_cmd->GetNewDoubleValue(newValue));
      G4cout << "Sample thickness changed to " << newValue << " ............................" << G4endl;
    }

  if(command == update_cmd)
    {
      myDetector->UpdateGeometry();
      G4cout << "Detector geometry updated ............................" << G4endl;
    }

}
