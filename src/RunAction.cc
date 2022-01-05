#include "RunAction.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "myAnalysisManager.hh"

#include "G4SDManager.hh"

RunAction::RunAction(){
	// Instantiate the custom analysis manager
	fAnalysisManager = new myAnalysisManager();
}

void RunAction::BeginOfRunAction(const G4Run*) {
	if (isMaster) G4Random::showEngineStatus();

	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->OpenFile();
	
	analysisManager->CreateNtuple("t0", "Event tree");
	analysisManager->CreateNtupleIColumn(0, "event");
	analysisManager->CreateNtupleIColumn(0, "det");
	analysisManager->CreateNtupleDColumn(0, "edep");
	analysisManager->CreateNtupleDColumn(0, "totsl");
	analysisManager->FinishNtuple(0);

	
	analysisManager->CreateNtuple("t1", "Detector tree");
	analysisManager->CreateNtupleIColumn(1, "event");
	analysisManager->CreateNtupleIColumn(1, "det");
	analysisManager->CreateNtupleIColumn(1, "strip");
	analysisManager->CreateNtupleIColumn(1, "layer");
	analysisManager->CreateNtupleDColumn(1, "edep");
	analysisManager->CreateNtupleIColumn(1, "pdg");
	analysisManager->FinishNtuple(1);
	

	analysisManager->CreateNtuple("t2", "Plane tree");
	analysisManager->CreateNtupleIColumn(2, "event");
	analysisManager->CreateNtupleDColumn(2, "x");
	analysisManager->CreateNtupleDColumn(2, "y");
	analysisManager->CreateNtupleDColumn(2, "z");
	analysisManager->CreateNtupleIColumn(2, "pdg");
	analysisManager->CreateNtupleDColumn(2, "ekin");
	analysisManager->CreateNtupleDColumn(2, "etot");
	analysisManager->CreateNtupleDColumn(2, "p");
	analysisManager->CreateNtupleDColumn(2, "px");//px -> componente x di p
	analysisManager->CreateNtupleDColumn(2, "py");
	analysisManager->CreateNtupleDColumn(2, "pz");
	analysisManager->FinishNtuple(2);
	

	analysisManager->CreateNtuple("t3", "Primary tree");
	analysisManager->CreateNtupleIColumn(3, "event");
	analysisManager->CreateNtupleDColumn(3, "x0");
	analysisManager->CreateNtupleDColumn(3, "y0");
	analysisManager->CreateNtupleDColumn(3, "z0");
	analysisManager->CreateNtupleDColumn(3, "px");
	analysisManager->CreateNtupleDColumn(3, "py");
	analysisManager->CreateNtupleDColumn(3, "pz");
	analysisManager->CreateNtupleDColumn(3, "ekin");
	analysisManager->CreateNtupleIColumn(3, "pdg");
	analysisManager->FinishNtuple(3);
	

	analysisManager->CreateNtuple("t4", "Strip tree");
	analysisManager->CreateNtupleIColumn(4, "event");
	analysisManager->CreateNtupleIColumn(4, "det");
	analysisManager->CreateNtupleIColumn(4, "strip");
	analysisManager->CreateNtupleDColumn(4, "edep");
	analysisManager->FinishNtuple(4);
	
	// Anton export tree data
	analysisManager->CreateNtuple("t5", "Anton export tree data");				// This Ntuple stores informations about charged particles
	//analysisManager->CreateNtupleIColumn(5, "det");								// Detector ID
	//analysisManager->CreateNtupleIColumn(5, "pdg");								// PDG code for the secondary part. type
	//analysisManager->CreateNtupleDColumn(5, "x");								// X position where charged secondary was recorded
	//analysisManager->CreateNtupleDColumn(5, "y");								// Y
	//analysisManager->CreateNtupleDColumn(5, "z");								// Z
	//analysisManager->CreateNtupleDColumn(5, "edep");							// Energy deposition
	analysisManager->FinishNtuple(5);
	
	// Pietro debug tree
	// This tree stores the energy deposition from pdgcode particle over the slen as a function of pos
	analysisManager->CreateNtuple("t6",	"Pietro debug tree");					// This Ntuple stores ...
	analysisManager->CreateNtupleIColumn(6, "event");
	analysisManager->CreateNtupleIColumn(6, "det");								// Detector ID
	analysisManager->CreateNtupleDColumn(6, "x");								// X
	analysisManager->CreateNtupleDColumn(6, "y");								// Y
	analysisManager->CreateNtupleDColumn(6, "z");								// Z
	analysisManager->CreateNtupleDColumn(6, "edep");							// Deposited energy
	analysisManager->CreateNtupleDColumn(6, "slen");							// Step length
	analysisManager->CreateNtupleIColumn(6, "pdg");								// PDG code
	analysisManager->FinishNtuple(6);
	
}

void RunAction::EndOfRunAction(const G4Run*) {
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	analysisManager->Write();
	analysisManager->CloseFile();

	if (isMaster) G4Random::showEngineStatus();
}