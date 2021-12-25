#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun) :
	G4UImessenger(),
	fAction(Gun)
{
	fGunDir = new G4UIdirectory("/profiler/gun/");
	fGunDir->SetGuidance("gun control");

	fPGunPosCmd = new G4UIcmdWith3VectorAndUnit("/profiler/gun/position", this);
	fPGunPosCmd->SetGuidance("Insert X Y and Z positions for the gamma beam distribution of the particle gun with respect to that of the \"World\"");
	fPGunPosCmd->SetParameterName("posx", "posy", "posz", false);
	fPGunPosCmd->SetDefaultUnit("cm");
	fPGunPosCmd->SetRange("posx>-25. && posx<25. && posy>-25. && posy<25. && posz>-25. && posz<25."); //sistemare il range rendendo il world autocalcolato // worldLength = 50 * cm;
	fPGunPosCmd->SetUnitCategory("Length");
	fPGunPosCmd->SetUnitCandidates("um mm cm");
	fPGunPosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fPGunPosZCmd = new G4UIcmdWithADoubleAndUnit("/profiler/gun/setZPosition", this);
	fPGunPosZCmd->SetGuidance("Set the z position of the particle gun");
	fPGunPosZCmd->SetParameterName("z", false);
	fPGunPosZCmd->SetUnitCategory("Length");
	fPGunPosZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fDefaultCmd = new G4UIcmdWithoutParameter("/profiler/gun/setDefault", this);
	fDefaultCmd->SetGuidance("set/reset kinematic defined in PrimaryGenerator");
	fDefaultCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fRndmCmd = new G4UIcmdWithADoubleAndUnit("/profiler/gun/gaussX", this);
	fRndmCmd->SetGuidance("Set the variance for the gaussian X lateral distribution of the beam");
	fRndmCmd->SetParameterName("sigma", false);
	fRndmCmd->SetRange("sigma>=0.");
	fRndmCmd->SetUnitCategory("Length");
	fRndmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

	fRndmCmdY = new G4UIcmdWithADoubleAndUnit("/profiler/gun/gaussY", this);
	fRndmCmdY->SetGuidance("Set the variance for the gaussian Y lateral distribution of the beam");
	fRndmCmdY->SetParameterName("sigma", false);
	fRndmCmdY->SetRange("sigma>=0.");
	fRndmCmdY->SetUnitCategory("Length");
	fRndmCmdY->AvailableForStates(G4State_PreInit, G4State_Idle);

	fImportBeamCmd = new G4UIcmdWithAString("/profiler/gun/importBeam", this);
	fImportBeamCmd->SetGuidance("Filename of the Kyle's ROOT file containing the beam");
	fImportBeamCmd->AvailableForStates(G4State_PreInit);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
	delete fDefaultCmd;
	delete fRndmCmd;
	delete fRndmCmdY;
	delete fGunDir;
	delete fPGunPosZCmd;
	delete fPGunPosCmd;
	delete fImportBeamCmd;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
	if (command == fDefaultCmd) fAction->SetDefaultKinematic();
	if (command == fRndmCmd)    fAction->SetGausSigmaX(fRndmCmd->GetNewDoubleValue(newValue));
	if (command == fRndmCmdY) fAction->SetGausSigmaY(fRndmCmdY->GetNewDoubleValue(newValue));
	if (command == fPGunPosZCmd) fAction->SetGunPosZ(fPGunPosZCmd->GetNewDoubleValue(newValue));
	if (command == fPGunPosCmd) fAction->SetPosition(fPGunPosCmd->GetNew3VectorValue(newValue));
	if (command == fImportBeamCmd) fAction->SetInputBeamFilename(newValue);
}
