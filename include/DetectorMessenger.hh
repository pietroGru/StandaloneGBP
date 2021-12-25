#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;

class DetectorMessenger: public G4UImessenger {
  
public:
  
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();
  
  virtual void SetNewValue(G4UIcommand*, G4String);
  
private:
  
  DetectorConstruction*      fDetector;
  G4UIdirectory*             fProfilerDir;
  G4UIdirectory*             fDetDir; 
  G4UIcmdWithAString*        fMatWorldCmd;
  G4UIcmdWithAString*        fMaterialCmd;
  G4UIcmdWithADoubleAndUnit* fThicknessCmd;    
  G4UIcmdWithADoubleAndUnit* fPitchCmd;    
  G4UIcmdWithADoubleAndUnit* fLengthCmd;    
  G4UIcmdWithAnInteger*      fNbStripsCmd;
  G4UIcmdWithAnInteger*      fNbLayersCmd;
  G4UIcmdWithADoubleAndUnit* fStripSpacingCmd;
  G4UIcmdWithADoubleAndUnit* fMetThicknessCmd;
};

#endif

