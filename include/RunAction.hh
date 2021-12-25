#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "myAnalysisManager.hh"
#include "G4AccumulableManager.hh"

class myAnalysisManager;

class RunAction : public G4UserRunAction {
public:
	RunAction();
	~RunAction() { delete fAnalysisManager; }

	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

private:
	myAnalysisManager* fAnalysisManager;
};

#endif

