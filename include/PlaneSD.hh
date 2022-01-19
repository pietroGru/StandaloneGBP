#ifndef PlaneSD_h
#define PlaneSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4VHitsCollection.hh"
#include "G4TouchableHistory.hh"
#include "planeHit.hh"

class G4HCofThisEvent;

class PlaneSD: public G4VSensitiveDetector {
public:
	PlaneSD(G4String name):
		G4VSensitiveDetector(name){
			collectionName.insert("inPartDownstreamCollection");
	};
	~PlaneSD(){};


	virtual void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	virtual void EndOfEvent(G4HCofThisEvent*);

private:
	planeHitsCollection* fHitsCollection;
	G4int fHCID;

	std::vector<G4int> recTrackID;
	std::vector<G4ThreeVector> recTrackPostStepPos;
};

#endif
