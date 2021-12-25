#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "mySD.hh"
#include "PlaneSD.hh"
//#include <iomanip>


DetectorConstruction::DetectorConstruction() :
	G4VUserDetectorConstruction(),
	fDefaultMaterial(0),
	fDetectorMessenger(0) {
	
	DefineMaterials();

	SetMatWorld("Galactic");
	SetMaterial("G4_Si");
	
	fMetalizationThickness = 100 * nm;
	fThickness = 100 * um;
	fPitch = 100 * um;
	fLength = 2 * cm;
	fNbOfStrips = 200;
	fNbOfLayers = 50;
	fSpacing = 10 * um;

	fDetectorMessenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction() { delete fDetectorMessenger; }

void DetectorConstruction::DefineMaterials() {
	G4NistManager* man = G4NistManager::Instance();

	// It seems that it is necessary to call the materials with this command for they to be available to the macro commands. (?)
	man->FindOrBuildMaterial("G4_Si");
	man->FindOrBuildMaterial("G4_AIR");
	
	G4Material* sapphire = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");

	G4double density = universe_mean_density;    //from PhysicalConstants.h
	G4double pressure = 3.e-18 * pascal;
	G4double temperature = 2.73 * kelvin;
	G4Material* Galactic = new G4Material("Galactic", 1., 1.008 * g / mole, density, kStateGas, temperature, pressure);

	fDefaultMaterial = Galactic;

	//Just a placeholder for now
	G4Element* Al = man->FindOrBuildElement("Al");
	fMetalizationMaterial = new G4Material("Aluminum", 13., 26.98 * g / mole, 2.700 * g / cm3);

	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


G4VPhysicalVolume* DetectorConstruction::ConstructVolumes() {
	// Cleanup old geometry
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	// TODO: optimize wordLenth to make run faster
	G4double worldLength = 50 * cm;
	// World
	G4Box* worldS = new G4Box("World", worldLength / 2., worldLength / 2., worldLength / 2.);
	G4LogicalVolume* worldL = new G4LogicalVolume(worldS, fMatWorld, "World");
	G4VPhysicalVolume* worldP = new G4PVPlacement(0, G4ThreeVector(), "World", worldL, 0, false, 0);

	// The detector block is for assumption stripped along the X axis.
	// Detector - sapphire
	G4Box* detectorS = new G4Box("Detector", fLength / 2., fLength / 2., fThickness / 2.);
	G4LogicalVolume* detectorL = new G4LogicalVolume(detectorS, fMaterial, "Detector");
	G4ThreeVector position(0, 0, 0);
	
	// Detector - metalization strip
	G4double mWidth = (fPitch - fSpacing);
	char tmp_mn[10]; sprintf(tmp_mn, "( %f um)", fThickness); G4String metName = tmp_mn;
	G4Box* metalizationS = new G4Box("Metalization", mWidth / 2., fLength / 2., fMetalizationThickness / 2.);
	G4LogicalVolume* metalizationL = new G4LogicalVolume(metalizationS, fMetalizationMaterial, "Metalization");
	//front strips
	G4ThreeVector strip_x0 =
		position																	// detector initial position
		+ G4ThreeVector(-(fNbOfStrips - 1) * fPitch / 2, 0, 0)						// minimum x range of the strip center
		+ G4ThreeVector(0, 0, (-fThickness + fMetalizationThickness) / 2);			// z offset to compensate for the thickness
	for (int i = 0; i < fNbOfStrips; i++) {
		new G4PVPlacement(0, strip_x0, metalizationL, "Strip", detectorL, false, i + 1, false);
		strip_x0 += G4ThreeVector(fPitch, 0., 0.);
	}

	// Detector - metalization plane (rear)
	G4Box* metalizationRearS = new G4Box("Metalization (rear)", fLength / 2., fLength / 2., fMetalizationThickness / 2.);
	G4LogicalVolume* metalizationRearL = new G4LogicalVolume(metalizationRearS, fMetalizationMaterial, "Metalization (rear)");
	G4ThreeVector strip_xN =
		position																									// detector initial position
		+ G4ThreeVector(0, 0, (fThickness - fMetalizationThickness) / 2);											// z offset to compensate for the thickness
	new G4PVPlacement(0, strip_xN, metalizationRearL, "Metalization (rear)", detectorL, false, 0, false);

	
	// Scoring plane
	//G4Box* scoringPlaneS = new G4Box("scoringPlaneS", worldLength / 2, worldLength / 2, fThickness / fNbOfLayers / 2);
	//G4LogicalVolume* scoringPlaneL = new G4LogicalVolume(scoringPlaneS, fDefaultMaterial, "scoringPlaneL");
	//G4VPhysicalVolume* scoringPlaneP = new G4PVPlacement(0, position + G4ThreeVector(0, 0, 2.5 * cm), "Scoring plane", scoringPlaneL, worldP, false, 2);
	//scoringPlaneL->SetVisAttributes(G4VisAttributes::GetInvisible());
	//G4cout << " pointer to part. scorer: " << scoringPlaneP << G4endl;
	// Scoring plane - charge scorer det0
	//G4Box* scoringPlaneChargeS = new G4Box("scoringPlaneChargeS", worldLength / 2, worldLength / 2, fThickness / fNbOfLayers / 2);
	//G4LogicalVolume* scoringPlaneChargeL = new G4LogicalVolume(scoringPlaneChargeS, fDefaultMaterial, "Charge scoring plane");
	//G4VPhysicalVolume* scoringPlaneChargeP = new G4PVPlacement(0, position + G4ThreeVector(0, 0, fThickness + (fThickness / fNbOfLayers / 2)), "Charge scoring plane", scoringPlaneChargeL, worldP, false, 2);
	//scoringPlaneChargeL->SetVisAttributes(G4VisAttributes::GetInvisible());
	//G4cout << " pointer to charge scorer: " << scoringPlaneChargeP << G4endl;

	// Detector placement
	// X-Axis detector (strips parallel to the x)
	G4VPhysicalVolume* detectorP_x = new G4PVPlacement(0, position, "Detector X", detectorL, worldP, false, 0);		
	// Y-Axis detector (strips parallel to the y)
	G4RotationMatrix* rot = new G4RotationMatrix(pi / 2, 0, 0);
	G4VPhysicalVolume* detectorP_y = new G4PVPlacement(rot, position + G4ThreeVector(0,0,2*cm), "Detector Y", detectorL, worldP, false, 1);
	
	// Set visualization attributes
	worldL->SetVisAttributes(G4VisAttributes::GetInvisible());

	G4VisAttributes* visAttr = new G4VisAttributes(G4Color(1, 0, 0, 0.5));
	detectorL->SetVisAttributes(visAttr);
	metalizationL->SetVisAttributes(G4Colour(0, 0, 1, 0.5));
	metalizationRearL->SetVisAttributes(G4Colour(0.2, 0, 0.8, 1.));	

	G4cout << " pointer to physiDet0: " << detectorP_x << G4endl;
	G4cout << " pointer to physiDet1: " << detectorP_y << G4endl;
	G4cout << " world material: " << worldL->GetMaterial()->GetName() << G4endl;
	return worldP;
}


G4VPhysicalVolume* DetectorConstruction::Construct() {
	return ConstructVolumes();
}

void DetectorConstruction::PrintParameters() {
	G4cout << "DetectorConstruction::PrintParameters"	<< G4endl
		<< "MatWorld:"									<< fMatWorld->GetName() << G4endl
		<< "Material:"									<< fMaterial->GetName() << G4endl
		<< "Metalization material:"						<< fMetalizationMaterial->GetName() << G4endl
		<< "Thickness:"									<< G4BestUnit(fThickness, "Length") << G4endl
		<< "Metalization thickness:"					<< G4BestUnit(fMetalizationThickness, "Length") << G4endl
		<< "Pitch:"										<< G4BestUnit(fPitch, "Length") << G4endl
		<< "Active strip lenght:\t"						<< G4BestUnit(fPitch - fSpacing, "Length") << G4endl
		<< "Length:"									<< G4BestUnit(fLength, "Length") << G4endl
		<< "Number of Strips:"							<< fNbOfStrips << G4endl
		<< "Number of Layers:"							<< fNbOfLayers << G4endl;
}

void DetectorConstruction::SetMatWorld(const G4String& material) {
	G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
	if (!pttoMaterial) {
		G4cout << " WARNING! Command refused. Material not found: " << material << G4endl;
		return;
	}
	fMatWorld = pttoMaterial;
	G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetMaterial(const G4String& material) {
	G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
	if (!pttoMaterial) {
		G4cout << " WARNING! Command refused. Material not found: " << material << G4endl;
		return;
	}
	fMaterial = pttoMaterial;
	G4RunManager::GetRunManager()->PhysicsHasBeenModified();
}

void DetectorConstruction::SetThickness(G4double val) {
	if (val <= DBL_MIN) {
		G4cout << " WARNING! Command refused. Thickness must be larger than " << DBL_MIN << G4endl;
		return;
	}
	fThickness = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetMetalizationThickness(G4double val) {
	if (val <= DBL_MIN || val >= fThickness) {
		G4cout	<< " WARNING! Command refused. Metalizationthickness must be larger than " << DBL_MIN
				<< " and lower than the detector thickness " << fThickness << G4endl;
		return;
	}
	fMetalizationThickness = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetPitch(G4double val) {
	if (val <= DBL_MIN) {
		G4cout << " WARNING! Command refused. Pitch must be larger than " << DBL_MIN << G4endl;
		return;
	}
	fPitch = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetStripSpacing(G4double val) {
	fSpacing = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetLength(G4double val) {
	if (val <= DBL_MIN) {
		G4cout << " WARNING! Command refused. Length must be larger than " << DBL_MIN << G4endl;
		return;
	}
	fLength = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetNbOfStrips(G4int val) {
	if (val < 1) {
		G4cout << " WARNING! Command refused. The number of strips must be larger than 0 " << G4endl;
		return;
	}
	fNbOfStrips = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetNbOfLayers(G4int val) {
	if (val < 1) {
		G4cout << " WARNING! Command refused. The number of layers must be larger than 0 " << G4endl;
		return;
	}
	fNbOfLayers = val;
	G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::ConstructSDandField() {
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	mySD* sensiDet = new mySD("sensiDet");					// Defines the sensitive GBP
	//sensiDet->SetVerboseLevel(2);
	sensiDet->SetNbOfStrips(fNbOfStrips);
	sensiDet->SetNbOfLayers(fNbOfLayers);
	sensiDet->SetThickness(fThickness);
	sensiDet->SetLength(fLength);
	//
	SDman->AddNewDetector(sensiDet);
	SetSensitiveDetector("Detector", sensiDet);

	//PlaneSD* planeSD = new PlaneSD("planeSD");			// Defines the sensitive scoring plane
	//SDman->AddNewDetector(planeSD);
	//SetSensitiveDetector("scoringPlaneL", planeSD);
	
	//PlaneSD* chargeScorer = new PlaneSD("chargeScorer");	// Defines the sensitive scoring plane #2 for charged
	//SDman->AddNewDetector(chargeScorer);
	//SetSensitiveDetector("Charge scoring plane", chargeScorer);
	
	G4cout << "detector costructed" << G4endl;
	if (fFieldMessenger.Get() == 0) {
		G4GlobalMagFieldMessenger* msg = new G4GlobalMagFieldMessenger(G4ThreeVector());
		//msg->SetVerboseLevel(1);
		G4AutoDelete::Register(msg);
		fFieldMessenger.Put(msg);
	}
}