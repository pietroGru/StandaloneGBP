#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys): 
  G4UImessenger(),
  fPhysicsList(pPhys) {

  fPhysDir = new G4UIdirectory("/profiler/phys/");
  fPhysDir->SetGuidance("physics list commands");

  fListCmd = new G4UIcmdWithAString("/profiler/phys/addPhysics", this);
  fListCmd->SetGuidance       ("Add modula physics list.");
  fListCmd->SetParameterName  ("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);
  fListCmd->SetToBeBroadcasted(false);

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/profiler/phys/stepMax", this);
  fStepMaxCmd->SetGuidance("Set max allowed step length");
  fStepMaxCmd->SetParameterName("mxStep", false);
  fStepMaxCmd->SetRange("mxStep>0.");
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PhysicsListMessenger::~PhysicsListMessenger() {
  delete fListCmd;
  delete fPhysDir;
  delete fStepMaxCmd;
}

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fListCmd) fPhysicsList->AddPhysicsList(newValue);
}
