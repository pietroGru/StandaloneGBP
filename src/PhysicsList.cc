#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "PhysListEmStandard.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ProcessManager.hh"

#include "StepMax.hh"


PhysicsList::PhysicsList() : G4VModularPhysicsList() {
	fMessenger = new PhysicsListMessenger(this);
	SetVerboseLevel(1);
	fEmPhysicsList = new PhysListEmStandard(fEmName = "local");
	G4LossTableManager::Instance();
	//G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(200 * eV, 100 * GeV);
	SetDefaultCutValue(1.0 * mm);
}

PhysicsList::~PhysicsList() {
	delete fMessenger;
	delete fEmPhysicsList;
}

void PhysicsList::ConstructParticle() {
	G4BosonConstructor      pBosonConstructor;      pBosonConstructor.ConstructParticle();
	G4LeptonConstructor     pLeptonConstructor;     pLeptonConstructor.ConstructParticle();
	G4MesonConstructor      pMesonConstructor;      pMesonConstructor.ConstructParticle();
	G4BaryonConstructor     pBaryonConstructor;     pBaryonConstructor.ConstructParticle();
	G4IonConstructor        pIonConstructor;        pIonConstructor.ConstructParticle();
	G4ShortLivedConstructor pShortLivedConstructor; pShortLivedConstructor.ConstructParticle();
}



void PhysicsList::ConstructProcess() {
	AddTransportation();
	fEmPhysicsList->ConstructProcess();
	G4EmParameters* param = G4EmParameters::Instance();
	param->SetBuildCSDARange(true);
	
	//AddDecay();
	//AddRadioactiveDecay();
	//AddStepMax();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
	// Step limitation seen as a process
	StepMax* stepMaxProcess = new StepMax(fMessenger);

	auto particleIterator = GetParticleIterator();
	particleIterator->reset();
	while ((*particleIterator)()) {
		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();

		if (stepMaxProcess->IsApplicable(*particle))
		{
			pmanager->AddDiscreteProcess(stepMaxProcess);
		}
	}
}


void PhysicsList::AddPhysicsList(const G4String& name) {
	if (verboseLevel > -1) G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
	
	if (name == fEmName) return;

	if (!(name == "local" ||
		name == "emstandard_opt0" ||
		name == "emstandard_opt1" ||
		name == "emstandard_opt2" ||
		name == "emstandard_opt3" ||
		name == "emstandard_opt4" ||
		name == "emstandardSS" ||
		name == "emstandardGS" ||
		name == "emstandardWVI" ||
		name == "emlowenergy" ||
		name == "emlivermore" ||
		name == "empenelope")) {
			G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << " is not defined" << G4endl;
			return;
	}

	fEmName = name;
	delete fEmPhysicsList;
	if (name == "local") {
		fEmPhysicsList = new PhysListEmStandard(name);
	}else if (name == "emstandard_opt0") {
		fEmPhysicsList = new G4EmStandardPhysics();
	}else if (name == "emstandard_opt1") {
		fEmPhysicsList = new G4EmStandardPhysics_option1();
	}else if (name == "emstandard_opt2") {
		fEmPhysicsList = new G4EmStandardPhysics_option2();
	}else if (name == "emstandard_opt3") {
		fEmPhysicsList = new G4EmStandardPhysics_option3();
	}else if (name == "emstandard_opt4") {
		fEmPhysicsList = new G4EmStandardPhysics_option4();
	}else if (name == "emstandardSS") {
		fEmPhysicsList = new G4EmStandardPhysicsSS();
	}else if (name == "emstandardGS") {
		fEmPhysicsList = new G4EmStandardPhysicsGS();
	}else if (name == "emstandardWVI") {
		fEmPhysicsList = new G4EmStandardPhysicsWVI();
	}else if (name == "emlowenergy") {
		fEmPhysicsList = new G4EmLowEPPhysics();
	}else if (name == "emlivermore") {
		fEmPhysicsList = new G4EmLivermorePhysics();
	}else if (name == "empenelope") {
		fEmPhysicsList = new G4EmPenelopePhysics();
	}else {
		G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << " is not defined" << G4endl;
	}
}

#include "G4Decay.hh"

void PhysicsList::AddDecay() {
	G4Decay* fDecayProcess = new G4Decay();

	auto particleIterator = GetParticleIterator();
	particleIterator->reset();

	while ((*particleIterator)()) {

		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();

		if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) { //@Pietro !particle->IsShortLived() come viene eseguito?

			pmanager->AddProcess(fDecayProcess);

			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager->SetProcessOrdering(fDecayProcess, idxPostStep);
			pmanager->SetProcessOrdering(fDecayProcess, idxAtRest);
		}
	}
}

#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"
#include "G4NuclideTable.hh"

void PhysicsList::AddRadioactiveDecay() {

	G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();

	radioactiveDecay->SetARM(true); // Atomic Rearangement

	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

	ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

	// mandatory for G4NuclideTable
	G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1 * picosecond);
}
