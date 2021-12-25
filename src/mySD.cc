#include "G4VProcess.hh"
#include "mySD.hh"
#include "myHit.hh"
#include "myAnalysisManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

void mySD::Initialize(G4HCofThisEvent* HCE) {
	// Create hits collection
	fHitsCollection = new myHitsCollection(GetName(), collectionName[0]);// GetName() -> SensitiveDetectorName
	
	// Cross-checks
	if (fHCID < 0) fHCID = GetCollectionID(0);
	if (fNbOfStrips >= STRIP_MAX) {	G4cout << " fNbOfStrips must be smaller than: " << STRIP_MAX << G4endl; exit(0); }
	if (fNbOfLayers >= LAYER_MAX) { G4cout << " fNbOfLayers must be smaller than: " << LAYER_MAX << G4endl; exit(0); }
	
	// Add this collection in hce
	HCE->AddHitsCollection(fHCID, fHitsCollection);
	
	// Verbose mode
	if (verboseLevel > 0) {
		//G4cout << " fNbOfStrips: " << fNbOfStrips << G4endl;
		//G4cout << " fNbOfLayers: " << fNbOfLayers << G4endl;
	}
}

G4bool mySD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
	G4double edep = aStep->GetTotalEnergyDeposit();						// Total energy deposition per step
	if (edep == 0.) return false;

	// Longitudinal profile of deposited energy (along z)
	G4StepPoint* prePoint = aStep->GetPreStepPoint();
	G4StepPoint* postPoint = aStep->GetPostStepPoint();
	G4ThreeVector P1 = prePoint->GetPosition();
	G4ThreeVector P2 = postPoint->GetPosition();
	G4ThreeVector point = P1 + G4UniformRand() * (P2 - P1);						// randomize point of energy deposition
	if (aStep->GetTrack()->GetDefinition()->GetPDGCharge() == 0.) point = P2;	// If the particle is neutral, the point of the step where it deposits energy is the last one

	G4double xshifted = point.x() + 0.5 * fLength;
	G4double yshifted = point.y() + 0.5 * fLength;
	G4double zshifted;

	auto touchable = aStep->GetPreStepPoint()->GetTouchableHandle();
	G4String physiName = touchable->GetVolume()->GetName();
	G4int detID;
	G4int stripID;
	G4int layerID;
	if (physiName == "Detector X"){
		detID = 0;
		zshifted = point.z() + 0.5 * fThickness;
		stripID = xshifted * (G4double)fNbOfStrips / fLength + 1;
	}else if (physiName == "Detector Y"){
		detID = 1;
		zshifted = point.z() - 2*CLHEP::cm + 0.5 * fThickness;
		stripID = yshifted * (G4double)fNbOfStrips / fLength + 1;
	}else{
		G4cout << " Unknown physical volume: " << physiName << G4endl;
		exit(0);
	}

	layerID = zshifted / (G4double)fThickness * (G4double)fNbOfLayers + 1;


	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4int pdgCode = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
	G4double sLeng = aStep->GetStepLength();
	// Dose evaluation
	analysisManager->FillNtupleIColumn(6, 0, detID);
	analysisManager->FillNtupleDColumn(6, 1, point.x());
	analysisManager->FillNtupleDColumn(6, 2, point.y());
	analysisManager->FillNtupleDColumn(6, 3, point.z());
	analysisManager->FillNtupleDColumn(6, 4, edep / CLHEP::keV);
	analysisManager->FillNtupleDColumn(6, 5, sLeng / CLHEP::um);
	analysisManager->FillNtupleIColumn(6, 6, pdgCode);
	analysisManager->AddNtupleRow(6);
	
	
	// Creates the hit object in the collection with the (det, strip, layer, edep) informations
	auto thisHit = new myHit(detID, stripID, layerID, pdgCode);
	thisHit->AddEdep(edep);
	thisHit->SetStepLength(sLeng);
	fHitsCollection->insert(thisHit);
	
	if (verboseLevel > 1) {
		thisHit->Print();
		//G4cout << "Z:" << zshifted << G4endl;
	}

	return true;
}

void mySD::EndOfEvent(G4HCofThisEvent*) {
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4RunManager* runManager = G4RunManager::GetRunManager();
	
	//Inizializza le variabili fEtot e fEdep (etot per strip) con 0
	// @Pietro. Perché usare auto quando ho un int?
	const G4int nbSDActive = 2;							// this variable is a shorthand to modify the sensitive detector number
	for (G4int det = 0; det < nbSDActive; det++) {
		fEtot[det] = 0;
		fSlentot[det] = 0;
		for (G4int strip = 0; strip < fNbOfStrips; strip++) {
			fEdep[det][strip] = 0;
		}
	}


	// fHitsCollection è un array (di un numero di entrate pari al numero di hit per evento) di 'array' (per la precisione di classi 'myHit') (di numero di entries pari al numero di variabili che ho definito io)
	size_t nHits = fHitsCollection->entries();
	G4int eventID = runManager->GetCurrentEvent()->GetEventID();

	if (verboseLevel > 1) {
		//G4cout << "Sensitive events:\t" << nHits << G4endl;
	}

	for (size_t hit = 0; hit < nHits; hit++) {
		auto myHit = (*fHitsCollection)[hit];

		G4double edep = myHit->GetEdep() / CLHEP::keV;
		if (!edep) continue; // se non c'è dep. di E nell'hit, allora passa al successivo

		if (verboseLevel > 1) myHit->Print();

		G4int detID = myHit->GetDetID();
		G4int stripID = myHit->GetStripID();
		G4int layerID = myHit->GetLayerID();
		G4int pdgCode = myHit->GetPDGCode();
		G4double stepLength = myHit->GetStepLength();
		
		if (stripID < 0 || layerID < 0) {
			char buffer[50];
			sprintf(buffer, "Strip/Layer ID not valid: %i, %i", stripID, layerID);
			G4Exception("mySD::EndOfEvent", "Negative array index", FatalException, buffer);
		}
		
		// Detector tree
		analysisManager->FillNtupleIColumn(1, 0, eventID);										// Event
		analysisManager->FillNtupleIColumn(1, 1, detID);										// Det
		analysisManager->FillNtupleIColumn(1, 2, stripID);										// Strip
		analysisManager->FillNtupleIColumn(1, 3, layerID);										// Layer
		analysisManager->FillNtupleDColumn(1, 4, edep);											// Deposited energy in this hit
		analysisManager->FillNtupleIColumn(1, 5, pdgCode);										// PDG code of the particle
		analysisManager->AddNtupleRow(1);
		
		fEtot[detID] += edep;
		fSlentot[detID] += powf(stepLength,3.);
		fEdep[detID][stripID] += edep;
	}

	// Event & Strip informations
	for (G4int det = 0; det < nbSDActive; det++) {
		// Record the energy deposition in the detector
		if (fEtot[det] == 0) { continue; }
		analysisManager->FillNtupleIColumn(0, 0, eventID);
		analysisManager->FillNtupleIColumn(0, 1, det);
		analysisManager->FillNtupleDColumn(0, 2, fEtot[det]);
		analysisManager->FillNtupleDColumn(0, 3, fSlentot[det]);
		analysisManager->AddNtupleRow(0);

		// Record the total energy deposition per strip
		for (G4int strip = 0; strip < fNbOfStrips; strip++) {
			if (fEdep[det][strip] == 0.) { continue; }
			// Strip tree
			analysisManager->FillNtupleIColumn(4, 0, eventID);
			analysisManager->FillNtupleIColumn(4, 1, det);
			analysisManager->FillNtupleIColumn(4, 2, strip);
			analysisManager->FillNtupleDColumn(4, 3, fEdep[det][strip]);
			analysisManager->AddNtupleRow(4);
		}
	}
}

