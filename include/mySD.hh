#ifndef mySD_h
#define mySD_h 1
#include "G4VHitsCollection.hh"
#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "myHit.hh"

#define STRIP_MAX 1000
#define LAYER_MAX 1000

class G4HCofThisEvent;

class mySD : public G4VSensitiveDetector {
public:
	mySD(G4String name):
		G4VSensitiveDetector(name),
		fHitsCollection(nullptr),
		fEdep(),
		fEtot(),
		fHCID(-1), fNbOfStrips(-1),
		fNbOfLayers(-1), fThickness(0),	fLength(0) {
		collectionName.insert("myCollection");
	};
	~mySD() {};
	
	inline void SetNbOfStrips(G4int val) { fNbOfStrips = val; }
	inline void SetNbOfLayers(G4int val) { fNbOfLayers = val; }
	inline void SetThickness(G4double val) { fThickness = val; }
	inline void SetLength(G4double val) { fLength = val; }
		
	virtual void Initialize(G4HCofThisEvent*);
	G4bool       ProcessHits(G4Step*, G4TouchableHistory*);
	virtual void EndOfEvent(G4HCofThisEvent*);

private:
	myHitsCollection* fHitsCollection;
	G4int fHCID;
	G4int fNbOfStrips;
	G4int fNbOfLayers;
	G4double fEtot[2];
	G4double fSlentot[2];
	G4double fEdep[2][STRIP_MAX]; // va definito dinamicamente ma come const usando fNbOfStrips -> così si spreca memoria

	G4double fThickness;
	G4double fLength;
};

#endif