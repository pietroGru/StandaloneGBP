#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector,
	PhysicsList* physics) :
	G4VUserActionInitialization(),
	fDetector(detector),
	fPhysics(physics)
{}

ActionInitialization::ActionInitialization(DetectorConstruction* detector) :
	G4VUserActionInitialization(),
	fDetector(detector) {
}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::BuildForMaster() const {
	SetUserAction(new RunAction());
}

void ActionInitialization::Build() const {
	PrimaryGeneratorAction* kinematics = new PrimaryGeneratorAction(fDetector);
	
	SetUserAction(kinematics);
	SetUserAction(new RunAction());

	EventAction* eventAction = new EventAction();

	SetUserAction(eventAction);
	SetUserAction(new TrackingAction());
	SetUserAction(new SteppingAction());
}

G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const {
	return new SteppingVerbose();
}
