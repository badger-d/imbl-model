#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;

class PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
  virtual ~PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:

  PrimaryGeneratorAction* Action;

  // Particle energy.
  G4UIcmdWithADoubleAndUnit* energy_cmd;

  // Isotropy.
  G4UIcmdWithAString* isotropy_cmd;

  // Particle origin.
  G4UIcmdWith3VectorAndUnit* origin_cmd;

};

#endif

