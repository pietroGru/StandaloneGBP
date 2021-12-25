#include "DetectorMessenger.hh"
#include <sstream>
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det) :
	G4UImessenger(),
	fDetector(Det) {
	fProfilerDir = new G4UIdirectory("/profiler/");
	fProfilerDir->SetGuidance(" detector control.");

	fDetDir = new G4UIdirectory("/profiler/det/");
	fDetDir->SetGuidance("detector construction commands");

	fMatWorldCmd = new G4UIcmdWithAString("/profiler/det/setMatWorld", this);
	fMatWorldCmd->SetGuidance("Set World Material name");
	fMatWorldCmd->SetParameterName("MatWorld", false);
	fMatWorldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fMatWorldCmd->SetToBeBroadcasted(false);

	fMaterialCmd = new G4UIcmdWithAString("/profiler/det/setMaterial", this);
	fMaterialCmd->SetGuidance("Set Material name");
	fMaterialCmd->SetParameterName("Material", false);
	fMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fMaterialCmd->SetToBeBroadcasted(false);

	fThicknessCmd = new G4UIcmdWithADoubleAndUnit("/profiler/det/setThickness", this);
	fThicknessCmd->SetGuidance("Set strip depth");
	fThicknessCmd->SetParameterName("Thickness", false);
	fThicknessCmd->SetRange("Thickness>0.");
	fThicknessCmd->SetUnitCategory("Length");
	fThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fThicknessCmd->SetToBeBroadcasted(false);

	fPitchCmd = new G4UIcmdWithADoubleAndUnit("/profiler/det/setPitch", this);
	fPitchCmd->SetGuidance("Set strip pitch");//@Pietro
	fPitchCmd->SetParameterName("Pitch", false);
	fPitchCmd->SetRange("Pitch>0.");
	fPitchCmd->SetUnitCategory("Length");
	fPitchCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fPitchCmd->SetToBeBroadcasted(false);

	fLengthCmd = new G4UIcmdWithADoubleAndUnit("/profiler/det/setLength", this);
	fLengthCmd->SetGuidance("Set strip length");
	fLengthCmd->SetParameterName("Length", false);
	fLengthCmd->SetRange("Length>0.");
	fLengthCmd->SetUnitCategory("Length");
	fLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fLengthCmd->SetToBeBroadcasted(false);

	fNbStripsCmd = new G4UIcmdWithAnInteger("/profiler/det/setNbOfStrips", this);
	fNbStripsCmd->SetGuidance("Set number of Strips.");
	fNbStripsCmd->SetParameterName("NbStrips", false);
	fNbStripsCmd->SetRange("NbStrips>0");
	fNbStripsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fNbStripsCmd->SetToBeBroadcasted(false);

	fNbLayersCmd = new G4UIcmdWithAnInteger("/profiler/det/setNbOfLayers", this);
	fNbLayersCmd->SetGuidance("Set number of Layers.");
	fNbLayersCmd->SetParameterName("NbLayers", false);
	fNbLayersCmd->SetRange("NbLayers>0");
	fNbLayersCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fNbLayersCmd->SetToBeBroadcasted(false);

	fStripSpacingCmd = new G4UIcmdWithADoubleAndUnit("/profiler/det/setStripSpacing", this);
	fStripSpacingCmd->SetGuidance("Set the spacing between a sensitive strip and another.");
	fStripSpacingCmd->SetParameterName("Spacing", false);
	fStripSpacingCmd->SetRange("Spacing>=0.");
	fStripSpacingCmd->SetUnitCategory("Length");
	fStripSpacingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fStripSpacingCmd->SetToBeBroadcasted(false);

	fMetThicknessCmd = new G4UIcmdWithADoubleAndUnit("/profiler/det/setMetThickness", this);
	fMetThicknessCmd->SetGuidance("Set the metalizzation thickness (front/rear have the same thickness)");
	fMetThicknessCmd->SetParameterName("Thickness", false);
	fMetThicknessCmd->SetRange("Thickness>=0.");
	fMetThicknessCmd->SetUnitCategory("Length");
	fMetThicknessCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
	fMetThicknessCmd->SetToBeBroadcasted(false);
}

DetectorMessenger::~DetectorMessenger() {
	delete fProfilerDir;
	delete fDetDir;
	delete fMatWorldCmd;
	delete fMaterialCmd;
	delete fThicknessCmd;
	delete fPitchCmd;
	delete fLengthCmd;
	delete fNbStripsCmd;
	delete fNbLayersCmd;
	delete fStripSpacingCmd;
	delete fMetThicknessCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
	if (command == fMatWorldCmd)  fDetector->SetMatWorld(newValue);
	if (command == fMaterialCmd)  fDetector->SetMaterial(newValue);
	if (command == fThicknessCmd) fDetector->SetThickness(fThicknessCmd->GetNewDoubleValue(newValue));
	if (command == fPitchCmd)     fDetector->SetPitch(fPitchCmd->GetNewDoubleValue(newValue));
	if (command == fLengthCmd)    fDetector->SetLength(fLengthCmd->GetNewDoubleValue(newValue));
	if (command == fNbStripsCmd)  fDetector->SetNbOfStrips(fNbStripsCmd->GetNewIntValue(newValue));
	if (command == fNbLayersCmd)  fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(newValue));
	if (command == fStripSpacingCmd)  fDetector->SetStripSpacing(fStripSpacingCmd->GetNewDoubleValue(newValue));
	if (command == fMetThicknessCmd)  fDetector->SetMetalizationThickness(fMetThicknessCmd->GetNewDoubleValue(newValue));
}
