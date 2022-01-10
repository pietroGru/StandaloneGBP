#include "planeHit.hh"
#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<planeHit>* planeHitAllocator;

planeHit::planeHit():
	G4VHit(),
	fEventID(-1),
	fTrackID(-1),
	fPDG(0),
	fEnergy(0.),
	fPosition(0),
	fMomentum(0.)
{}

planeHit::planeHit(const planeHit& right):
	G4VHit(),
	fEventID(right.fEventID),
	fTrackID(right.fTrackID),
	fPDG(right.fPDG),
	fEnergy(right.fEnergy),
	fPosition(right.fPosition),
	fMomentum(right.fMomentum)
{}

planeHit::planeHit(G4int trackID, G4int pdg, G4double ekin, G4ThreeVector position, G4ThreeVector momentum) :
	G4VHit(),
	fTrackID(trackID),
	fPDG(pdg),
	fEnergy(ekin),
	fPosition(position),
	fMomentum(momentum)
{}

planeHit::planeHit(G4int eventID, G4int trackID, G4int pdg, G4double ekin, G4ThreeVector position, G4ThreeVector momentum):
	G4VHit(),
	fEventID(eventID),
	fTrackID(trackID),
	fPDG(pdg),
	fEnergy(ekin),
	fPosition(position),
	fMomentum(momentum)
{}

planeHit::~planeHit() {}

const planeHit& planeHit::operator=(const planeHit& right) {
	fEventID = right.fEventID;
	fTrackID = right.fTrackID;
	fPDG = right.fPDG;
	fEnergy = right.fEnergy;
	fPosition = right.fPosition;
	fMomentum = right.fMomentum;
	return *this;
}

G4bool planeHit::operator<(const planeHit& right) {
	std::cout << "(" << fTrackID << "<" << right.fTrackID << ")\t";
	return fTrackID > right.fTrackID;
}


//G4bool planeHit::operator==(const planeHit& right) const {
//	return (fEventID == right.fEventID && fTrackID == right.fTrackID);
//}

void planeHit::Print() {
	G4cout << " EventID: " << fEventID
		<< " TrackID: " << fTrackID
		<< " PDG: " << fPDG
		<< " Energy: " << fEnergy
		<< " Position: " << fPosition
		<< " Momentum: " << fMomentum
		<< G4endl;
}