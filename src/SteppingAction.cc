#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "Randomize.hh"

void SteppingAction::UserSteppingAction(const G4Step*) {}