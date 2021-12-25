#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "Randomize.hh"

class G4Event;
class DetectorConstruction;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

	PrimaryGeneratorAction(DetectorConstruction*);
	~PrimaryGeneratorAction();

	//Gaussian distributed beam. Functions
	void SetGausSigmaX(G4double val) { fRndmBeam = val; beamMode = 1; }
	void SetGausSigmaY(G4double val) { fRndmBeamY = val; beamMode = 1; }
	void SetDefaultKinematic();

	void SetInputBeamFilename(G4String name) { ImportBeamFromFile(name); beamMode = 2; }

	virtual void SetPosition(G4double x0, G4double y0, G4double z0);
	void SetPosition(G4ThreeVector pos) { pGunPosition = pos; }
	void SetGunPosZ(G4double val) { pGunPosition.setZ(val); }
	
	virtual void GeneratePrimaries(G4Event*);

	G4ParticleGun* GetParticleGun() { return fParticleGun; }

private:
	G4ParticleGun* fParticleGun;
	DetectorConstruction* fDetector;
	G4ThreeVector			   pGunPosition;
	PrimaryGeneratorMessenger* fGunMessenger;
	
	G4int beamMode;																// This variable sets the beam type.
	// 0-point source with fixed pos/mom;
	// 1-gaussian 2D beam;
	// 2-realistic beam imported from root file;

	// Gaussian distributed beam. Variables
	G4double                   fRndmBeam;
	G4double                   fRndmBeamY;
	// Functions
	G4ThreeVector GaussianBeam_pos(G4ThreeVector center);


	// Import beam distribution from file
	// Variables
	std::vector<G4ThreeVector> fGunKyle_pos, fGunKyle_mom;
	std::vector<G4double> fGunKyle_energy, fGunKyle_weight;
	std::vector<G4int> fGunKyle_pdg;
	int tracksEntries;
	// Functions
	void ImportBeamFromFile(G4String beamFilename); //import kyle file as a whole
	G4ThreeVector ImportedBeamFile_pos(int eventNb);
	G4ThreeVector ImportedBeamFile_mom(int eventNb);
	G4double ImportedBeamFile_energy(int eventNb);
	G4double ImportedBeamFile_weight(int eventNb);
	G4int ImportedBeamFile_pdg(int eventNb);


	// Import beam distribution using probability density function (PDF)
	// Variables
	//std::vector<G4RandGeneral> fGunRandDistr;
	//std::vector<G4double> fGunRandDistr_range;
	//std::vector<G4double> fGunRandDistr_mean;
	//std::vector<G4double> fGunRandDistr_x0;
	// Functions
	//void ImportBeamDistPDF();
	//G4ThreeVector ImportedBeam_pos();
	//G4ThreeVector ImportedBeam_mom();
};

#endif

