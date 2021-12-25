//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm1/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Track.hh"
#include "StackingAction.hh"
#include "myAnalysisManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
	: G4UserStackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* track)
{
	
	//keep primary particle
	if (track->GetParentID() == 0) return fUrgent;
	
	//
	//energy spectrum of secondaries
	//
	G4double energy = track->GetKineticEnergy();
	G4double charge = track->GetDefinition()->GetPDGCharge();
	G4ThreeVector position = track->GetPosition();
	
	
	auto touchable = track->GetTouchableHandle();
	G4String physiName = touchable->GetVolume()->GetName();
	G4int detID=-1;
	if (physiName == "World") {
		return fUrgent;
	}else if (physiName == "Detector X" || physiName == "Strip" || physiName == "Metalization (rear)") {
		detID = 0;
	}
	else if (physiName == "Detector Y" || physiName == "Strip" || physiName == "Metalization (rear)") {
		detID = 1;
	}
	else{
		G4cout << " Unknown physical volume: " << physiName << G4endl;
		return fUrgent; // exit(0);
	}
	
	if (detID != -1 && charge != 0.){
		G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
		analysisManager->FillNtupleIColumn(6, 0, detID);
		analysisManager->FillNtupleDColumn(6, 1, position.x());
		analysisManager->FillNtupleDColumn(6, 2, position.y());
		analysisManager->FillNtupleDColumn(6, 3, position.z());
		analysisManager->FillNtupleDColumn(6, 4, energy / CLHEP::keV);
		analysisManager->AddNtupleRow(6);
	}
	
	return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
