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

  ion_cham_cmd = new G4UIcmdWithAString("/geometry/ion_cham",this);
  ion_cham_cmd->SetGuidance("Selects inclusion of insulated box");
  ion_cham_cmd->SetGuidance("Choice : on(default), off");
  ion_cham_cmd->SetParameterName("choice",true);
  ion_cham_cmd->SetDefaultValue("on");
  ion_cham_cmd->SetCandidates("on off");
  ion_cham_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sample_cmd = new G4UIcmdWithAString("/geometry/sample",this);
  sample_cmd->SetGuidance("Selects inclusion of sample");
  sample_cmd->SetGuidance("Choice : on(default), off");
  sample_cmd->SetParameterName("choice",true);
  sample_cmd->SetDefaultValue("on");
  sample_cmd->SetCandidates("on off");
  sample_cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sample_pos_cmd = new G4UIcmdWith3VectorAndUnit("/geometry/sample_pos", this);
  sample_pos_cmd->SetGuidance("Set position of the sample.");
  sample_pos_cmd->SetParameterName("ox","oy","oz",true,true);
  sample_pos_cmd->SetDefaultUnit("mm");

  update_cmd = new G4UIcmdWithoutParameter("/geometry/update",this);
  update_cmd->SetGuidance("Update detector geometry.");
  update_cmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  update_cmd->SetGuidance("if you changed geometrical value(s).");
  update_cmd->AvailableForStates(G4State_Idle);

}


DetectorMessenger::~DetectorMessenger()
{

  delete ion_cham_cmd;
  delete sample_cmd;
  delete sample_pos_cmd;
  delete update_cmd;

}

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{


  if(command == ion_cham_cmd)
    {
      myDetector->Set_Ion_Cham_Flag(newValue);
      G4cout << "### Ion chamber included." << newValue << G4endl;
    }

  if(command == sample_cmd)
    {
	  myDetector->Set_Sample_Flag(newValue);
	  G4cout << "### Sample included." << newValue << G4endl;
    }

  if(command == sample_pos_cmd)
    {
      myDetector->Set_Sample_Pos(sample_pos_cmd->GetNew3VectorValue(newValue));
      G4cout << "Sample position changed to " << newValue << " ............................" << G4endl;
    }

  if(command == update_cmd)
    {
      myDetector->UpdateGeometry();
      G4cout << "Detector geometry updated ............................" << G4endl;
    }

}
