#ifndef planeHit_h
#define planeHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class planeHit:public G4VHit {
public:
	planeHit();
	planeHit(const planeHit& right);
	planeHit(G4int trackID, G4int pdg, G4double ekin, G4ThreeVector position, G4ThreeVector momentum);
	planeHit(G4int eventID, G4int trackID, G4int pdg, G4double ekin, G4ThreeVector position, G4ThreeVector momentum);
	~planeHit(); //perché virtual? Non esiste il decostruttore di planeHit!

	const planeHit& operator=(const planeHit& right);
	G4bool operator<(const planeHit& right);
	//G4bool operator==(const planeHit& right) const;

	inline void* operator new(size_t);
	inline void operator delete(void* aHit);

	virtual void Print();

	inline void SetEventID(G4int id) { fEventID = id; }
	inline void SetTrackID(G4int id) { fTrackID = id; }
	inline void SetPDG(G4int pdg) { fPDG = pdg; }
	inline void SetEnergy(G4double ekin) { fEnergy = ekin; }
	inline void SetPosition(G4ThreeVector pos) { fPosition = pos; }
	inline void SetMomentumDirection(G4ThreeVector mom) { fMomentum = mom; }

	G4int GetEventID() { return fEventID; }
	G4int GetTrackID() { return fTrackID; }
	G4int GetPDG() { return fPDG; }
	G4double GetEnergy() { return fEnergy; }
	G4ThreeVector GetPosition() { return fPosition; }
	G4ThreeVector GetMomentumDirection() { return fMomentum; }

private:
	G4int fEventID, fTrackID;
	G4int fPDG;
	G4ThreeVector fPosition, fMomentum;
	G4double fEnergy;
};


using planeHitsCollection = G4THitsCollection<planeHit>;
extern G4ThreadLocal G4Allocator<planeHit>* planeHitAllocator;

inline void* planeHit::operator new(size_t) {
	if (!planeHitAllocator) { planeHitAllocator = new G4Allocator<planeHit>; }
	return (void*)planeHitAllocator->MallocSingle();
}

inline void planeHit::operator delete(void* aHit) {
	planeHitAllocator->FreeSingle((planeHit*)aHit);
}

#endif
