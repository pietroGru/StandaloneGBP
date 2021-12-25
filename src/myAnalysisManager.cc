#include "myAnalysisManager.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"

void myAnalysisManager::Book() {
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetFileName(fFileName);
	analysisManager->SetVerboseLevel(1);
	analysisManager->SetActivation(true);
}
/*
void myAnalysisManager::AddNtuple(G4int id){
	if (id > maxrecordNTupleNb) {
		char buffer[100];
		sprintf(buffer, "The tuple with ID %i exceeded the maximum allower number of recordables.", id);
		G4Exception("myAnalysisManager::AddNtuple", "New tuple instance over max limit", FatalException, buffer);
	}
	else {
		tupleIDStore[nbTuple] = id;
		nbTuple++;
	}
}
*/