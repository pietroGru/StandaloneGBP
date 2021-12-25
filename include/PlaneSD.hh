#ifndef PlaneSD_h
#define PlaneSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"

class PlaneSD: public G4VSensitiveDetector {
  
public:
  
  PlaneSD(G4String name): G4VSensitiveDetector(name) {};
  ~PlaneSD() {};
  
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
};

#endif
