#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "PlaneSD.hh"
#include "planeHit.hh"
#include "myAnalysisManager.hh"


void PlaneSD::Initialize(G4HCofThisEvent* HCE){
	// Create hits collection
	fHitsCollection = new planeHitsCollection(GetName(), collectionName[0]);
	// Cross-checks
	if (fHCID < 0) fHCID = GetCollectionID(0);
	// Add this collection in hce
	HCE->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool PlaneSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
	// Get the post step point and determine the kinetic energy of the particle there
	G4StepPoint* postPoint = aStep->GetPostStepPoint();
	// As well as the position of this last point
	G4ThreeVector pPointPos = postPoint->GetPosition();
	G4ThreeVector pPointMom = postPoint->GetMomentumDirection();

	// Get the track ID associated with this step
	G4Track* aTrack = aStep->GetTrack();
	G4int trackID = aTrack->GetTrackID();

	
	G4int pdg = aTrack->GetParticleDefinition()->GetPDGEncoding();
	G4double ekin = postPoint->GetKineticEnergy();

	// Particle is entering the detector
	if (pPointMom.z() > 0) {
		auto thisHit = new planeHit(trackID, pdg, ekin, pPointPos, pPointMom);
		fHitsCollection->insert(thisHit);

		if (verboseLevel > 1) {
			thisHit->Print();
		}
	}
	// This part can be used later
	//aTrack->SetTrackStatus(fStopAndKill);
	//aTrack->SetTrackStatus(fKillTrackAndSecondaries);
	return true;
}

G4bool PlaneSD_compareFcn(planeHit* left, planeHit* right) {
	return left->GetTrackID() < right->GetTrackID();
}


void PlaneSD::EndOfEvent(G4HCofThisEvent* HCE){
	G4RunManager* runManager = G4RunManager::GetRunManager();
	G4int eventID = runManager->GetCurrentEvent()->GetEventID();
	// Sort the collection depending upon the trackID number
	std::sort(fHitsCollection->GetVector()->begin(), fHitsCollection->GetVector()->end(), PlaneSD_compareFcn);

	// Remove duplicates from hits one is no interested in
	auto vector = *fHitsCollection->GetVector();
	G4int trackID_bak=1;
	for (size_t i = 1; i < vector.size(); i++) {
		auto hit = vector[i];
		G4int trackID = hit->GetTrackID();
		G4double postStepPosZ = hit->GetPosition().z();

		if (trackID == trackID_bak) {
			// Debug. Prints the whole vector
			//auto vectorDebug = vector;
			//G4cout << "Original size was: " << vectorDebug.size() << G4endl;
			//for (size_t ii = 0; ii < vectorDebug.size(); ii++) {
			//	auto hitDebug = vectorDebug[ii];
			//	G4int trackIDDebug = hitDebug->GetTrackID();
			//	G4cout << trackIDDebug << ", ";
			//}
			//G4cout << G4endl;
			//----processing---
			auto hitPrev = vector[i-1];
			G4double postStepPosZPrev = hitPrev->GetPosition().z();
			if (postStepPosZ <= postStepPosZPrev) {
				vector.erase(vector.begin() + i);
			}
			else {
				vector.erase(vector.begin() + i-1);
			}
			i--;
			//----after removal----
			// Debug. Prints the whole vector
			//vectorDebug = vector;
			//G4cout << "Actual size is: " << vectorDebug.size() << G4endl;
			//for (size_t ii = 0; ii < vectorDebug.size(); ii++) {
			//	auto hitDebug = vectorDebug[ii];
			//	G4int trackIDDebug = hitDebug->GetTrackID();
			//	G4cout << trackIDDebug << ", ";
			//}
			//G4cout << G4endl;
		}
		trackID_bak = trackID;
	}

	// Save info in the file
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	size_t nHits = vector.size();
	for (size_t i = 0; i < nHits; i++) {
		auto hit = vector[i];
		analysisManager->FillNtupleIColumn(2, 0, eventID);
		analysisManager->FillNtupleIColumn(2, 1, hit->GetTrackID());
		analysisManager->FillNtupleIColumn(2, 2, hit->GetPDG());
		analysisManager->FillNtupleDColumn(2, 3, hit->GetPosition().x() / CLHEP::mm);
		analysisManager->FillNtupleDColumn(2, 4, hit->GetPosition().y() / CLHEP::mm);
		analysisManager->FillNtupleDColumn(2, 5, hit->GetPosition().z() / CLHEP::mm);
		analysisManager->FillNtupleDColumn(2, 6, hit->GetMomentumDirection().x());
		analysisManager->FillNtupleDColumn(2, 7, hit->GetMomentumDirection().y());
		analysisManager->FillNtupleDColumn(2, 8, hit->GetMomentumDirection().z());
		analysisManager->FillNtupleDColumn(2, 9, hit->GetEnergy() / CLHEP::GeV);
		analysisManager->AddNtupleRow(2);
	}
}