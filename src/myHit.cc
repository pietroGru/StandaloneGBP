#include "myHit.hh"
#include "DetectorConstruction.hh"
#include "G4Box.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<myHit>* myHitAllocator;

myHit::myHit() :
	G4VHit(),
	fDetID(-1),
	fStripID(-1),
	fLayerID(-1),
	fEdep(0.),
	fParticleType(-1),
	fStepLength(0.)
{}

myHit::myHit(G4int detID, G4int stripID, G4int layerID) :
	G4VHit(),
	fDetID(detID),
	fStripID(stripID),
	fLayerID(layerID),
	fEdep(0.),
	fParticleType(-1),
	fStepLength(0.)
{}

myHit::myHit(G4int detID, G4int stripID, G4int layerID, G4double stepLength) :
	G4VHit(),
	fDetID(detID),
	fStripID(stripID),
	fLayerID(layerID),
	fEdep(0.),
	fParticleType(-1),
	fStepLength(stepLength)
{}

myHit::myHit(G4int detID, G4int stripID, G4int layerID, G4int particlePDGCode) :
	G4VHit(),
	fDetID(detID),
	fStripID(stripID),
	fLayerID(layerID),
	fEdep(0.),
	fParticleType(particlePDGCode)
{}

myHit::myHit(G4int detID, G4int stripID, G4int layerID, G4int particlePDGCode, G4double stepLength) :
	G4VHit(),
	fDetID(detID),
	fStripID(stripID),
	fLayerID(layerID),
	fEdep(0.),
	fParticleType(particlePDGCode),
	fStepLength(stepLength)
{}

myHit::~myHit() {}

myHit::myHit(const myHit& right) :
	G4VHit(),
	fDetID(right.fDetID),
	fStripID(right.fStripID),
	fLayerID(right.fLayerID),
	fEdep(right.fEdep),
	fParticleType(right.fParticleType)
{}

const myHit& myHit::operator=(const myHit& right) {
	fDetID = right.fDetID;
	fStripID = right.fStripID;
	fLayerID = right.fLayerID;
	fEdep = right.fEdep;
	fParticleType = right.fParticleType;
	return *this;
}

G4bool myHit::operator==(const myHit& right) const {
	return (fDetID == right.fDetID && fStripID == right.fStripID && fLayerID == right.fLayerID);
}

void myHit::Print() {
	G4cout << " DetID: " << fDetID
		<< " StripID: " << fStripID
		<< " LayerID: " << fLayerID
		<< " PDG: " << fParticleType
		<< " Edep[keV]: " << fEdep / keV
		<< " Step len.[um]: " << fStepLength / um
		<< G4endl;
}
