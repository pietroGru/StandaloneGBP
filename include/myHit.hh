#ifndef myHit_h
#define myHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class myHit : public G4VHit{
public:
	myHit();
	myHit(G4int detID, G4int stripID, G4int layerID);
	myHit(G4int detID, G4int stripID, G4int layerID, G4double stepLength);
	myHit(G4int detID, G4int stripID, G4int layerID, G4int particlePDGCode);
	myHit(G4int detID, G4int stripID, G4int layerID, G4int particlePDGCode, G4double stepLength);
	myHit(const myHit& right);
	~myHit(); //perché virtual? Non esiste il decostruttore di myHit!

	const myHit& operator=(const myHit& right);
	G4bool      operator==(const myHit& right) const;

	inline void* operator new(size_t);
	inline void operator delete(void* aHit);

	virtual void Print();

	inline void SetDetID(G4int    id) { fDetID = id; }
	inline void SetStripID(G4int    id) { fStripID = id; }
	inline void SetLayerID(G4int    id) { fLayerID = id; }
	inline void SetEdep(G4double de) { fEdep = de; }
	inline void SetStepLength(G4double sl) { fStepLength = sl; }
	inline void SetParticleType(G4int pdg) { fParticleType = pdg; }
	inline void AddEdep(G4double de) { fEdep += de; }
	inline void AddStepLength(G4double sl) { fStepLength += sl; }
	G4int    GetDetID()   const { return fDetID; }
	G4int    GetStripID() const { return fStripID; }
	G4int    GetLayerID() const { return fLayerID; }
	G4int    GetPDGCode() const { return fParticleType; }
	G4double GetEdep()    const { return fEdep; }
	G4double GetStepLength()    const { return fStepLength; }
	//G4int GetChargeCouples() const { return fChargeCouples; }

private:
	G4int fDetID;
	G4int fStripID;
	G4int fLayerID;
	G4double fEdep;
	G4int fParticleType;
	G4double fStepLength;
	//Gint fChargeCouples;
};


using myHitsCollection = G4THitsCollection<myHit>;
extern G4ThreadLocal G4Allocator<myHit>* myHitAllocator;

inline void* myHit::operator new(size_t) {
	if (!myHitAllocator) { myHitAllocator = new G4Allocator<myHit>; }
	return (void*)myHitAllocator->MallocSingle();
}

inline void myHit::operator delete(void* aHit) {
	myHitAllocator->FreeSingle((myHit*)aHit);
}

#endif
