#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;

class PrimaryGeneratorMessenger : public G4UImessenger {
public:
	PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
	~PrimaryGeneratorMessenger();
	virtual void SetNewValue(G4UIcommand*, G4String);
private:
	PrimaryGeneratorAction* fAction;
	G4UIdirectory* fGunDir;
	G4UIcmdWithoutParameter* fDefaultCmd;
	G4UIcmdWithADoubleAndUnit* fRndmCmd;
	G4UIcmdWithADoubleAndUnit* fRndmCmdY;
	G4UIcmdWithADoubleAndUnit* fPGunPosZCmd;
	G4UIcmdWith3VectorAndUnit* fPGunPosCmd;
	G4UIcmdWithAString* fImportBeamCmd;
};

#endif

