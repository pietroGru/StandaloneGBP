#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "PlaneSD.hh"
#include "myAnalysisManager.hh"

G4bool PlaneSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
	G4Track* aTrack = aStep->GetTrack();

	//aTrack->SetTrackStatus(fStopAndKill);
	//aTrack->SetTrackStatus(fKillTrackAndSecondaries);
	G4int trackID = aTrack->GetTrackID();
	if (trackID == 1) return true; //a che serve sta cosa?

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	G4RunManager* runManager = G4RunManager::GetRunManager();

	const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();

	analysisManager->FillNtupleIColumn(2, 0, runManager->GetCurrentEvent()->GetEventID());
	analysisManager->FillNtupleDColumn(2, 1, aTrack->GetPosition().x());
	analysisManager->FillNtupleDColumn(2, 2, aTrack->GetPosition().y());
	analysisManager->FillNtupleDColumn(2, 3, aTrack->GetPosition().z());
	analysisManager->FillNtupleIColumn(2, 4, aParticle->GetPDGcode());
	analysisManager->FillNtupleDColumn(2, 5, aParticle->GetKineticEnergy());
	analysisManager->FillNtupleDColumn(2, 6, aParticle->GetTotalEnergy());
	analysisManager->FillNtupleDColumn(2, 7, aParticle->GetTotalMomentum());
	analysisManager->FillNtupleDColumn(2, 8, aParticle->GetMomentumDirection().x());
	analysisManager->FillNtupleDColumn(2, 9, aParticle->GetMomentumDirection().y());
	analysisManager->FillNtupleDColumn(2, 10, aParticle->GetMomentumDirection().z());
	analysisManager->AddNtupleRow(2);
	
	return true;
}
