#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {

public:

	DetectorConstruction();
	~DetectorConstruction();

	void SetMatWorld(const G4String&);
	void SetMaterial(const G4String&);
	void SetThickness(G4double);
	void SetMetalizationThickness(G4double);
	void SetPitch(G4double);
	void SetStripSpacing(G4double);
	void SetLength(G4double);
	void SetNbOfStrips(G4int);
	void SetNbOfLayers(G4int);

	virtual G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructVolumes();
	virtual void ConstructSDandField();

	G4Material* GetMatWorld() { return fMatWorld; };
	G4Material* GetMaterial() { return fMaterial; };
	G4double    GetThickness() { return fThickness; };
	G4double    GetMetalizationThickness() { return fMetalizationThickness; };
	G4double    GetPitch() { return fPitch; };
	G4double    GetStripSpacing() { return fSpacing; };
	G4double    GetLength() { return fLength; };
	G4int       GetNbOfStrips() { return fNbOfStrips; }
	G4int       GetNbOfLayers() { return fNbOfLayers; }

	void PrintParameters();

private:
	G4double    fLength;							// Detector length
	G4double    fThickness;							// Detector thickness
	G4double    fMetalizationThickness;				// Thickness of the metal layer (Al)
	G4double    fPitch;								// Distance between a strip and anohter
	G4double    fSpacing;							// Spacing between two strips with metalization
	G4int       fNbOfStrips;						// Number of strips

	G4Material* fMatWorld;							// World material
	G4Material* fMaterial;							// Material for the detector
	G4Material* fMetalizationMaterial;				// Material for the metalized layer

	G4int       fNbOfLayers;						// Number of layers (unused)
	G4Material* fDefaultMaterial;					// Default material

	DetectorMessenger* fDetectorMessenger;
	G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
private:
	void DefineMaterials();
};

#endif
