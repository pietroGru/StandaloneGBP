#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class PhysicsListMessenger: public G4UImessenger {
  
public:
  
  PhysicsListMessenger(PhysicsList* p = 0);
  ~PhysicsListMessenger();
  
  virtual void SetNewValue(G4UIcommand*, G4String);
  inline G4double GetMaxChargedStep() const { return fMaxChargedStep; }

private:
  
  PhysicsList* fPhysicsList;
  
  G4UIdirectory*      fPhysDir;        
  G4UIcmdWithAString* fListCmd;
  G4UIcmdWithADoubleAndUnit* fStepMaxCmd;
  G4double fMaxChargedStep;
};

#endif

