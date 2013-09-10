#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction*);
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  DetectorConstruction* myDetector;
  G4UIdirectory*        geometryDir;

  // Commands to switch components on/off.
  G4UIcmdWithAString* detector_cmd;
  G4UIcmdWithAString* sample_cmd;

  // Command to update the sample thickness.
  G4UIcmdWithADoubleAndUnit* sampleX_cmd;
  
  // Command to update the geometry.
  G4UIcmdWithoutParameter* update_cmd;

};

#endif

