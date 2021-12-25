#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

G4int debugger;

int main(int argc, char** argv) {
	//detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = nullptr;
	if (argc == 1) ui = new G4UIExecutive(argc, argv);

	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4VSteppingVerbose::SetInstance(new SteppingVerbose);
#ifdef G4MULTITHREADED
	G4MTRunManager* runManager = new G4MTRunManager;
	G4int nThreads = std::min(G4Threading::G4GetNumberOfCores(), 4);
	if (argc == 3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
	runManager->SetNumberOfThreads(nThreads);
	G4cout << "===== profiler is started with " << runManager->GetNumberOfThreads() << " threads =====" << G4endl;
#else
	G4RunManager* runManager = new G4RunManager;
#endif

	DetectorConstruction* det = new DetectorConstruction;
	PhysicsList* phys = new PhysicsList;

	runManager->SetUserInitialization(det);
	runManager->SetUserInitialization(phys);
	runManager->SetUserInitialization(new ActionInitialization(det, phys));

	G4VisManager* visManager = nullptr;
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	if (ui) { // interactive mode
		visManager = new G4VisExecutive;
		visManager->Initialize();
		UImanager->ApplyCommand("/control/execute vis.mac");
		ui->SessionStart();
		delete ui;
	}
	else { // batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command + fileName);
	}

	// @Pietro - for not to close the terminal in batch mode
	do {
		G4cout << '\n' << "Press the Enter key to close." << G4endl;
	} while (G4cin.get() != '\n');

	delete visManager;
	delete runManager;
}
