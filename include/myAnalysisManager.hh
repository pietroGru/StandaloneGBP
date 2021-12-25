#ifndef myAnalysisManager_h
#define myAnalysisManager_h 1

#include "globals.hh"
#include "g4root.hh" //L'inclusione di sto header dice che l'analisi viene fatta usando root

class myAnalysisManager {
public:

  myAnalysisManager(): fFileName("run") { Book(); }
  ~myAnalysisManager() { delete G4AnalysisManager::Instance(); }
private:
  
  void Book();
  G4String fFileName;
};

#endif

