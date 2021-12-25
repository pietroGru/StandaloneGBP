#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "myAnalysisManager.hh"
#include "g4root.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det) :
	G4VUserPrimaryGeneratorAction(),
	fDetector(det),
	fRndmBeam(0.),
	pGunPosition(0),
	beamMode(-1){
		fParticleGun = new G4ParticleGun(1);
		fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(22)); // Photon
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
		fParticleGun->SetParticleEnergy(5. * GeV);

		fGunMessenger = new PrimaryGeneratorMessenger(this);
		SetDefaultKinematic();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
	delete fParticleGun;
	delete fGunMessenger;
}

void PrimaryGeneratorAction::SetDefaultKinematic() {
	fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(22)); // Photon
	fParticleGun->SetParticleEnergy(5. * GeV);
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
	fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -0.5 * fDetector->GetThickness()));
}

void PrimaryGeneratorAction::SetPosition(G4double x0, G4double y0, G4double z0){
	G4ThreeVector ppGun_pos = G4ThreeVector(x0, y0, z0);
	fParticleGun->SetParticlePosition(ppGun_pos);
	pGunPosition = ppGun_pos;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
	// The z position of the particle gun is moved so that the distance from the first surface of the detector is 1 nm.
	// This is just for convenience, because I wonder if precision errors are accumulating because of the identical z-positions
	// of both the detector and the particle source.
	auto pGunM = fParticleGun->GetParticleMomentumDirection();
	G4ThreeVector ppGun_pos = pGunPosition;
	G4ThreeVector ppGun_mom = G4ThreeVector(pGunM.getX(), pGunM.getY(), pGunM.getZ());
	G4double ppGun_energy = fParticleGun->GetParticleEnergy();

	switch (beamMode){
		case 0:
			//pointlike source
			break;

		case 1:
			// Gaussian distributed beam
			ppGun_pos = GaussianBeam_pos(pGunPosition);
			fParticleGun->SetParticlePosition(ppGun_pos);
			break;

		case 2:
			// Beam from file
			G4int evId = anEvent->GetEventID();
			if (evId > tracksEntries) {
				// events simulated are bigger than the size of the ptarmigan
				fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
				fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1.));
			}
			ppGun_pos += ImportedBeamFile_pos(evId) * mm;
			ppGun_mom = ImportedBeamFile_mom(evId) * GeV; /* /c */
			ppGun_energy = ImportedBeamFile_energy(evId) * GeV;
			G4int ppGun_pdg = ImportedBeamFile_pdg(evId);
			G4int ppGun_weight = ImportedBeamFile_weight(evId);
			//
			fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(ppGun_pdg));
			fParticleGun->SetNumberOfParticles(ppGun_weight);
			fParticleGun->SetParticleEnergy(ppGun_energy);
			fParticleGun->SetParticlePosition(ppGun_pos);
			fParticleGun->SetParticleMomentumDirection(ppGun_mom.unit());
			break;
	}

	fParticleGun->GeneratePrimaryVertex(anEvent);

	// Record information about primary particle
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->FillNtupleIColumn(3, 0, anEvent->GetEventID());
	analysisManager->FillNtupleDColumn(3, 1, ppGun_pos.getX() / mm);
	analysisManager->FillNtupleDColumn(3, 2, ppGun_pos.getY() / mm);
	analysisManager->FillNtupleDColumn(3, 3, ppGun_pos.getZ() / mm);
	analysisManager->FillNtupleDColumn(3, 4, ppGun_mom.getX() / GeV); // GeV/c
	analysisManager->FillNtupleDColumn(3, 5, ppGun_mom.getY() / GeV); // GeV/c
	analysisManager->FillNtupleDColumn(3, 6, ppGun_mom.getZ() / GeV); // GeV/c
	analysisManager->FillNtupleDColumn(3, 7, fParticleGun->GetParticleEnergy() / GeV);
	analysisManager->FillNtupleIColumn(3, 8, fParticleGun->GetParticleDefinition()->GetPDGEncoding());
	analysisManager->AddNtupleRow(3);
}

// Gaussian beam mode
// Return the three-vector for the particle gun position distributed with sigmaX, Y over the center mean value
G4ThreeVector PrimaryGeneratorAction::GaussianBeam_pos(G4ThreeVector center) {
	G4ThreeVector ppGun_pos = center;
	if (fRndmBeam > 0. && fRndmBeamY > 0.) { // randomize the beam along XY with different variances
		//@Pietro
		if (fRndmBeam > fDetector->GetLength()) {
			fRndmBeam = fDetector->GetLength();
		}
		else if (fRndmBeamY > 2. * cm) {
			fRndmBeamY = 2. * cm;
		}
		ppGun_pos += G4ThreeVector(G4RandGauss::shoot(0., fRndmBeam), G4RandGauss::shoot(0., fRndmBeamY), 0.);
	}
	else if (fRndmBeam > 0.) { // randomize the beam along X with fRndmBeam variance
		if (fRndmBeam > fDetector->GetLength()) fRndmBeam = fDetector->GetLength();
		ppGun_pos += G4ThreeVector(G4RandGauss::shoot(0., fRndmBeam), 0., 0.);
	}
	else if (fRndmBeamY > 0.) { // randomize the beam along Y with fRndmBeamY variance
		if (fRndmBeamY > 2. * cm) fRndmBeamY = 2. * cm;
		ppGun_pos += G4ThreeVector(0., G4RandGauss::shoot(0., fRndmBeamY), 0.);
	}
	return ppGun_pos;
}



// Input beam mode 1 - File
// Read the ROOT file with input beam data and fill local variables with the info.s contained in the Tracks tree
void PrimaryGeneratorAction::ImportBeamFromFile(G4String beamFilename) {
	// Create (or get) analysis reader
	G4AnalysisReader* analysisReader = G4AnalysisReader::Instance();
	analysisReader->SetVerboseLevel(0);
	if(analysisReader->GetVerboseLevel() > 0)	G4cout << "Request to import beam spec.s from the file: " << beamFilename << G4endl;
	
	analysisReader->SetFileName(beamFilename);
	G4String hName[] = { "x", "y", "z", "px", "py", "pz", "energy", "weight", "pdg" };
	const int varNb = std::size(hName);
	G4double vars[varNb-1];
	G4int varPDG;
	G4int trackEntry = 0;
	
	G4int ntupleId = analysisReader->GetNtuple("Tracks");
	if (ntupleId >= 0) {
		G4cout << "Setting branch addresses...";
		bool check = false;
		for (int i = 0; i < varNb-1; i++) {
			check |= analysisReader->SetNtupleDColumn(hName[i], vars[i]);
			G4cout << hName[i] << ", ";
		}
		check |= analysisReader->SetNtupleIColumn(hName[varNb-1], varPDG);
		G4cout << hName[varNb - 1] << "...";
		
		if (!check) {
			G4Exception("PrimaryGeneratorAction::ImportBeamDist", "Unable to set branch addresses", FatalException, "");
		}
		else {
			G4cout << "passed!" << G4endl;
		}

		// Read ntuple
		if (analysisReader->GetVerboseLevel() > 0) G4cout << "Reading ntuple...";
		while (analysisReader->GetNtupleRow()) {
			tracksEntries++;
			if (analysisReader->GetVerboseLevel() > 1) {
				G4cout << trackEntry << "th entry: ";
				for (int i = 0; i < (varNb-1); i++) {
					G4cout << hName[i] << ": " << vars[i] << ", ";
				}
				G4cout << hName[varNb-1] << ": " << varPDG << G4endl;
			}
			fGunKyle_pos.push_back(G4ThreeVector(vars[0], vars[1], vars[2]));
			fGunKyle_mom.push_back(G4ThreeVector(vars[3], vars[4], vars[5]));
			fGunKyle_energy.push_back(vars[6]);
			fGunKyle_weight.push_back(vars[7]);
			fGunKyle_pdg.push_back(varPDG);
		}
		if (analysisReader->GetVerboseLevel() > 0) G4cout << "finished! " << trackEntry << " entries." << G4endl;
	}
	else {
		G4Exception("PrimaryGeneratorAction::ImportBeamDistFile", "Unable to load TTree", FatalException, "Tracks");
	}

	// Set beam mode to "import from file"
	beamMode = 1;
}
// Return the position vector for the imported particle at TTree entry number #eventNb
G4ThreeVector PrimaryGeneratorAction::ImportedBeamFile_pos(int eventNb) {
	return fGunKyle_pos[eventNb];
}
// Return the momentum vector for the imported particle at TTree entry number #eventNb
G4ThreeVector PrimaryGeneratorAction::ImportedBeamFile_mom(int eventNb) {
	return fGunKyle_mom[eventNb];
}
// Return the energy for the imported particle at TTree entry number #eventNb
G4double PrimaryGeneratorAction::ImportedBeamFile_energy(int eventNb) {
	return fGunKyle_energy[eventNb];
}
// Return the MC weight for the imported particle at TTree entry number #eventNb
G4double PrimaryGeneratorAction::ImportedBeamFile_weight(int eventNb) {
	return fGunKyle_weight[eventNb];
}
// Return the pdg encoding for the imported particle at TTree entry number #eventNb
G4int PrimaryGeneratorAction::ImportedBeamFile_pdg(int eventNb) {
	return fGunKyle_pdg[eventNb];
}








/*
// TODO
// Input beam mode 2 - PDF
void PrimaryGeneratorAction::ImportBeamDistPDF() {
}

G4ThreeVector PrimaryGeneratorAction::ImportedBeam_pos() {
	G4double tpos[3] = {};
	for (int i = 0; i < 3; i++) {
		if (fGunRandDistr_range[i] != 0) {
			tpos[i] = fGunRandDistr_x0[i] + fGunRandDistr[i].shoot() * fGunRandDistr_range[i];
		}
		else {
			tpos[i] = fGunRandDistr_mean[i];
		}
	}
	return G4ThreeVector(tpos[0], tpos[1], tpos[2]);
}

G4ThreeVector PrimaryGeneratorAction::ImportedBeam_mom() {
	G4double tmom[3] = {};
	for (int i = 3; i < 6; i++) {
		if (fGunRandDistr_range[i] != 0) {
			tmom[i-3] = fGunRandDistr_x0[i] + fGunRandDistr[i].shoot() * fGunRandDistr_range[i];
		}
		else {
			tmom[i-3] = fGunRandDistr_mean[i];
		}
	}
	return G4ThreeVector(tmom[0], tmom[1], tmom[2]);
}
*/