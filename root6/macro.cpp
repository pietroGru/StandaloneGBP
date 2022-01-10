/*
ROOT Analysis macro for the StandaloneGBP
Tree description
t0 - Event tree
    0 - event           It stores the event identifier.
    1 - det             It is either 0 or 1. It can be used to discriminate between the two sensitive detectors.
    2 - edep            The energy deposition (keV) in the hit.
t1 - Detector tree
It is like a digitization tree. The ditigization is made by discretizing the hit positions over the SD elements (X-strips/Z-layers).
    0 - event           It stores the event identifier.
    1 - det             It is either 0 or 1. It can be used to discriminate between the two sensitive detectors.
    2 - strip           It stores the strip number [0,fNbOfStrips] where the hit occurred.
    3 - layer           It stores the layer number [0,fNbOfLayers] where the hit occurred.
    4 - edep            The energy deposition (keV) in the hit.

t2 - Plane tree
The scoring plane is optionally placed after the two detectors. It records informations about the particles arriving at the place. 
	0 - event           It stores the event identifier.
    1 - x               The x position of the hit.
    2 - y               The y position of the hit.
    3 - z               The z position of the hit.
    4 - pdg             The PDG code for the particle type hitting the plane.
    5 - ekin            The kinetic energy of the particle hitting the plane.
    6 - p               Momentum.
    7 - px              Momentum x.
    8 - py              Momentum y.
    9 - pz              Momentum z.

t3 - Primary tree
This tree contains info about the primary action, that it the initial gamma beam. 
	0 - event           It stores the event identifier.
    1 - x0              The x initial position of the particle gun.
    2 - y0              The y initial position of the particle gun.
    3 - z0              The z initial position of the particle gun.

t4 - Strip tree
This tree contains the total deposited energy per strip per detector.
	0 - event           It stores the event identifier.
    1 - det             It is either 0 or 1. It can be used to discriminate between the two sensitive detectors.
    2 - strip           It stores the strip number [0,fNbOfStrips] where the hit occurred.
    4 - edep            The cumulative energy deposition (keV) in the strip (layer-independent).

t5 - Charge secondary tracks
This tree stores the event when a charged particle (of secondary type) hits the sensitive detector.
    0 - det             It is either 0 or 1. It can be used to discriminate between the two sensitive detectors.
    1 - x               The position.
    2 - y               The position.
    3 - z               The position.
    4 - edep            The energy deposition (keV) in the hit.

t6 - Debug tree - dose evaluation
This tree stores the energy deposition from pdgcode particle over the slen as a function of pos
    0 - det             It is either 0 or 1. It can be used to discriminate between the two sensitive detectors.
    1 - x               The position.
    2 - y               The position.
    3 - z               The position.
    4 - edep            Energy deposition during the step
    5 - slen            Step length
    6 - pdg             PDG encoding for the particle type
*/
// Counters for progress indicator and minimal impact on routine
#include "tools.cpp"
int verbosityLevel = 1;                                                                         // Verbosity levels: 0-none, 1-info, 2-debug

//Detector geometry
double gausBSigma[2] = {375.0, 375.0}; //um
double detThickness = 100.0; //um 


TString ifilename = "kyle/kyleBeam_9_50_0.root";
TFile* input = new TFile(ifilename, "read");
TTree* eventTree        = (TTree*)input->Get("t0");
TTree* detTree          = (TTree*)input->Get("t1");
TTree* planeTree        = (TTree*)input->Get("t2");
TTree* primaryTree      = (TTree*)input->Get("t3");
TTree* stripTree        = (TTree*)input->Get("t4");
TTree* chargeTree       = (TTree*)input->Get("t5");
TTree* debugTree        = (TTree*)input->Get("t6");

// Utilities
    void loadAnotherFile(TString filename);                                                     // Load a new file
    TString whichDetVerbose(int detNb);                                                         // Returns either 'upstream' or 'downstream' depending if i=0 or i=1
    TString whichDetShort(int detNb);                                                           // Returns 'det0' or 'det1'
    TString whichConfigIs();                                                                    // Returns the configuration prefix for file
    int getConfig(TString);
    inline int getConfig(){ return getConfig(ifilename); }
// Primary
    std::vector<TH1D> beamMonitor(int bins = 400, float range = 2, bool autobin = false);       // Primary particle monitor
    std::vector<TH2D> primaryProfile(int pdg=11);                                               // Number XY profile of the particles. Filter for particles with given pdg code
    std::vector<TH1D> primarySpectrum();                                                        // Energy spectrum of primaries
    std::vector<TH1D> primarySpectrum(int pdg, double range = 10.0);                            // Energy spectrum of primaries with the given pdg code. Return 2-array with 0-pdgfiltered, 1-all part.s
    std::vector<TH1D> primarySpectrumLogLog(int pdg, double lowRange = 1E-9,                    // Energy spectrum of primaries with the given pdg code, with variable size bins. Return 2-array with 0-pdgfiltered, 1-all part.s
        double highRange = 16.0, int binsNb = 100);
    std::vector<TH2D> primaryEnergyProfile();                                                   // Energy XY profile of the particles. Each point (x,y) is weightned with the corresponding energy
// Energy deposition plots
    std::vector<TH1D> edepSpectrum(int bins, float emin, float emax);                           // Energy spectrum of depositions for the upstream/downstream detector
    std::vector<TH1D> edepSpectrumPartitionedByPrimaryEnergies(int bins=1000,                   // Energy spectrum of depositions for the upstream/downstream detector in a window Ea-Eb, partitioned by primary energies. Energy in keV
    double emin = 0.0 /*GeV*/, double emax= 6.0 /*GeV*/, double baseRatio = -1E6 /*keV*/);
    std::vector<TH1D> edepSpectrumChg(int bins, double a, double b);                            // Energy spectrum of depositions for selected particles. Return a 2-array with the spectrum of energy deposited in upstream/downstream detector in the range [a,b] with #bins. Sel is a selector for the values \xi={"5.0", "7.0", "10.0"}
    std::vector<TH1D> edepSpectrumChg();                                                        // Energy spectrum of depositions for selected particles. Default energy range 0-1000 and 1000 bins.
    std::vector<TH1D> eSpectrumStrip(int stripNb = 100, float range = 5.);                      // Energy spectrum of depositions in a given strip number for some det#
    double totalEdep(int detNb=0);                                                              // Energy deposited in detNb in keV
    double* eStrip_array(int stripNb = 100);                                                    // 2-array of energy deposited in the strip 100 for upstream/downstream detectors in keV
    double eStrip(int detNb=0, int stripNb = 100);                                              // Energy deposited in the strip 'stripNb' of detector 'detNb' in keV
    // double eStrip(int detNb=0, int stripNb = 100, float range = 1.);                         // Energy spectrum of depositions in a given strip number. Return energy deposited in keV
    std::vector<TH1D> eLongitudinal_array(int nbBins = 100);                                    // Longitudinal energy deposition profile. Custom binning available. Return 2-array with det0 and det1 plots
    TH1D eLongitudinal(int detNb);                                                              // Longitudinal energy deposition profile for detNb
    TH1D eLongitudinal(int detNb, int nbBins);                                                  // Longitudinal energy deposition profile for detNb with custom binning.
    std::vector<TH1D> edepProject(double meshSize = 0.100);                                     // Energy deposition in the transverse plane. Return 4-array [0]det0x, [1]det0y, [2]det1x, [3]ddet1y
    std::vector<TH1D> edepProject(bool filter, double meshSize = 0.100);                        // Energy deposition in the transverse plane for primaries with a selection filter (i.e. charge/neutral). Return 4-array [0]det0x, [1]det0y, [2]det1x, [3]ddet1y
    std::vector<TH1D> edepProjectStripTree();                                                   // Energy deposition profile in the strips. It uses the strip tree data. Return a 2-array for det0 and 1
    std::vector<TH1D> edepProjectStripTree(double pE0, double pE1);                             // Energy deposition profile in the strips given primaries in selected energy range [pE0, pE1] GeV. Return a 2-array for det0 and 1
    std::vector<TH1D> edepStripChg();                                                           // Energy deposition profile in the strips with comparison with deposition from charged primaries, for the primaries in the energy range [enLow, enHigh]. It uses the strip tree data. Return a 2-array for det0 and 1
    std::vector<TH1D> edepProjectPrimEnPartitioned(double enLow=0 /*GeV*/,                      // Energy depositions in the strips separated between the contributions with different energies
        double enHigh=10.0 /*GeV*/, double baseRatio = 10.0);
    void plot_edepProjectPrimEnPartitioned(double enLow=0 /*GeV*/,                              //TEMPDESCRIPTION-Plot function edepProjectPrimEnPartitioned
        double enHigh=10.0 /*GeV*/, double baseRatio = 10.0, 
        TString path = "plots/",TString format = ".pdf");
// Dose
    std::vector<TH2D> doseXY(double meshSize = 0.100 /*mm*/);                                   // Energy & dose 2D map using a mesh of size LxLx100 um3. Return a 4-array of TH2D with [0]emapdet0, [1]emapdet1, [2]dmapdet0, [3]dmapdet1
    std::vector<TH2D> doseXY(bool neutral, double meshSize = 0.100 /*mm*/);                     // Energy & dose 2D map using a mesh of size LxLx100 um3. Contributions from only photons to the dose. Return a 4-array of TH2D with [0]emapdet0, [1]emapdet1, [2]dmapdet0, [3]dmapdet1
// Profile reconstruction    
    std::vector<TH1I> fitHits(double threshold, int* fitRange);                                 // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Threshold is in keV. Return a 2-array with fit for upstream and downstream det.s
    std::vector<TH1I> fitHits();                                                                // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Overloading
    double* fitHitsChi2(double threshold, int* fitRange);                                       // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Threshold is in keV, fitRange is a 2-array with stripL, stripH where the fit is ranged.
    double* fitHitsChi2(double threshold = 1.);                                                 // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Default threshold of 1keV and range 90,110.
    std::vector<TGraph2D> fitVSThrRng();                                                        // 2D-graph of the chisquare as a function of threshold and range
    void plot_fitVSThrRng();                                                                    //TEMPDESCRIPTION-Plot fitVSThrRng()
// Other methods + unsorted
    void overlayReconstructedBeamProfile();                                                     // Superimpose the beam profile of the primaries with the reconstructed one
    int createReport(TString path = "plots/", TString format = ".pdf",                          // Save plots in the folder plots
        TString ofilename = "out.root");
// Other
    void toolsCumulativeProbabilityXYGaussian(double events = 1E5, float sig = 0.375,           // This function explains the choose for the radius used in the average dose calculation.
        float mult = 1.);



void loadAnotherFile(TString name){
    TFile* input = new TFile(name, "read");
    eventTree        = (TTree*)input->Get("t0");
    detTree          = (TTree*)input->Get("t1");
    planeTree        = (TTree*)input->Get("t2");
    primaryTree      = (TTree*)input->Get("t3");
    stripTree        = (TTree*)input->Get("t4");
    chargeTree       = (TTree*)input->Get("t5");
    debugTree        = (TTree*)input->Get("t6");
    ifilename = name;
}
// Returns either 'upstream' or 'downstream' depending if i=0 or i=1
TString whichDetVerbose(int i){
    TString detName;
    if(i==0){
        detName = "upstream";
    }else if(i==1){
        detName = "downstream";
    }else{
        detName = "???";
    }
    return detName;
}
// Returns 'det0' or 'det1'
TString whichDetShort(int i){
    return TString::Format("det%i", i);
}
// Returns the configuration prefix for file
TString whichConfigIs(){
    TString spec;
    // Number of photons
    if(ifilename.Contains("_7_")){
        spec = "7_";
    }else if(ifilename.Contains("_9_")){
        spec = "9_";
    }
    // Sigmaconfig
    if(ifilename.Contains("_r12_")){
        spec += "r12_";
    }else if(ifilename.Contains("_r13_")){
        spec += "r13_";
    }
    return spec;
}
// Log in a base base
double log(double base, double x){
    return std::log(x) / std::log(base);
}



// Primary particle monitor
std::vector<TH1D> beamMonitor(int bins = 400, float range = 2, bool autobin = false){
    TH1D result[3];
    result[0] = TH1D("beamX", "Beam monitor;x [mm];counts", bins, -range, range);
    result[1] = TH1D("beamY", "Beam monitor;y [mm];counts", bins, -range, range);
    result[2] = TH1D("beamZ", "Beam monitor;z [mm];counts", 110, -110, 0);
    primaryTree->Project("beamX", "x0");
    primaryTree->Project("beamY", "y0");
    primaryTree->Project("beamZ", "z0");
    if(verbosityLevel > 0){
        TCanvas* beamMonitorCanvas = new TCanvas("beamMonitorCanvas", "Primary particle monitor", 1600, 400);
        beamMonitorCanvas->Divide(3);
        for(int i=0; i<3; i++){
            beamMonitorCanvas->cd(i+1);
            result[i].DrawClone();
        }
    }
    std::vector<TH1D> retVar;
    for(int i=0; i<3; i++){
        retVar.push_back(result[i]);
    }
    return retVar;
}
// Number XY profile of the particles.
std::vector<TH2D> primaryProfile(int pdg=11){
    TH2D* primaryHist = new TH2D("primaryHist", "Transverse primaries profile / BX; x [mm]; y [mm]; counts", 200, -10., 10., 200, -10., 10.);
    TH2D* primaryChgHist = new TH2D("primaryChgHist", "Transverse primaries profile for |pdg|=11 / BX; x [mm]; y [mm]; counts", 200, -10., 10., 200, -10., 10.);
    primaryTree->Project("primaryHist", "y0:x0");
    primaryTree->Project("primaryChgHist", "y0:x0", TString::Format("pdg==%i || pdg==%i", pdg, -pdg));
    primaryHist->SetDrawOption("LEGO1");
    if(verbosityLevel > 0){
        TCanvas* primaryProfileCanvas = new TCanvas("primaryProfileCanvas", "Number XY profile of the particles");
        primaryHist->DrawClone("LEGO1");
        new TCanvas();
        primaryChgHist->DrawClone("surf1");
    }

    std::vector<TH2D> retVar;
    retVar.push_back(*primaryHist);
    retVar.push_back(*primaryChgHist);
    delete primaryHist, primaryChgHist;
    return retVar;
}
// Energy spectrum of primaries
std::vector<TH1D> primarySpectrum(){
    TH1D* primaryeSpectrHist = new TH1D("primaryeSpectrHist", "Energy spectrum of primaries; energy [GeV]; counts", 100, 0, 10);
    primaryTree->Project("primaryeSpectrHist", "ekin");
    if(verbosityLevel > 0){
        TCanvas* primarySpectrumCanvas = new TCanvas("primarySpectrumCanvas", "Energy spectrum of primaries");
        primarySpectrumCanvas->SetGridx(); primarySpectrumCanvas->SetGridy();
        primarySpectrumCanvas->SetLogx(); primarySpectrumCanvas->SetLogy();
        primaryeSpectrHist->DrawClone();
    }
    
    std::vector<TH1D> retVar;
    retVar.push_back(*primaryeSpectrHist);
    delete primaryeSpectrHist;
    return retVar;
}
// Energy spectrum of primaries with the given pdg code. Return 2-array with 0-pdgfiltered, 1-all part.s
std::vector<TH1D> primarySpectrum(int pdg, double range = 10.0){
    // cout << "breakpoint A. pars: " << pdg << '\t' << range << endl;
    TH1D* primaryeSpectrHist = new TH1D("primaryeSpectrHist", "Energy spectrum of primaries; energy [GeV]; counts", range*10.0, 0, range);
    TH1D* primaryeSpectrPDGHist = new TH1D("primaryeSpectrPDGHist", TString::Format("Energy spectrum of primaries with pdg == %i || pdg== %i; energy [GeV]; counts", pdg, -pdg), range*10.0, 0, range);
    // cout << "breakpoint A2. pars: " <<  range*10.0 << '\t' <<  0 << '\t' << range << endl;
    primaryTree->Project("primaryeSpectrHist", "ekin");
    // cout << "breakpoint A22. pars: " <<  range*10.0 << '\t' <<  0 << '\t' << range << endl;
    primaryTree->Project("primaryeSpectrPDGHist", "ekin", TString::Format("pdg==%i || pdg==%i",pdg, -pdg));
    // cout << "breakpoint A3" << endl;
    primaryeSpectrPDGHist->SetLineColor(kRed);
    if(verbosityLevel>0){
        TCanvas* primarySpectrumPDGCanvas = new TCanvas("primarySpectrumPDGCanvas", TString::Format("Energy spectrum of pdg==%i primaries", pdg));
        primarySpectrumPDGCanvas->SetGridx(); primarySpectrumPDGCanvas->SetGridy();
        primarySpectrumPDGCanvas->SetLogx(); primarySpectrumPDGCanvas->SetLogy();
        primaryeSpectrHist->DrawClone();
        primaryeSpectrPDGHist->DrawClone("same");
        cout << TString::Format("Total entries: %i, Chg. entries: %i\n", primaryeSpectrHist->GetEntries(), primaryeSpectrPDGHist->GetEntries());
    }
    // cout << "breakpoint B" << endl;
    std::vector<TH1D> retVar;
    retVar.push_back(*primaryeSpectrPDGHist);
    retVar.push_back(*primaryeSpectrHist);
    // cout << "breakpoint C" << endl;
    return retVar;
}
// Energy spectrum of primaries with the given pdg code, with variable size bins. Return 2-array with 0-pdgfiltered, 1-all part.s
std::vector<TH1D> primarySpectrumLogLog(int pdg, double lowRange = 1E-9, double highRange = 16.0, int binsNb = 100){
    // Creating the bins array 
    double* bins = (double*)malloc((binsNb+1)*sizeof(double));
    bins[0] = 0;
    double N = (binsNb-1)/( 1 - (TMath::Log(lowRange)/TMath::Log(highRange)));
    for(int i=0; i<binsNb; i++){
        bins[binsNb-i] = pow(highRange, 1 - i/N);
        // cout << bins[binsNb-i] << '\t';
    }
    // cout << "breakpoint A. pars: " << pdg << '\t' << range << endl;
    TH1D* primaryeSpectrHist = new TH1D("primaryeSpectrLogBinsHist", "Energy spectrum of primaries; energy [GeV]; counts", binsNb, bins);
    TH1D* primaryeSpectrPDGHist = new TH1D("primaryeSpectrPDGLogBinsHist", TString::Format("Energy spectrum of primaries with pdg == %i || pdg== %i; energy [GeV]; counts", pdg, -pdg), binsNb, bins);
    // cout << "breakpoint A2. pars: " <<  range*10.0 << '\t' <<  0 << '\t' << range << endl;
    primaryTree->Project("primaryeSpectrLogBinsHist", "ekin");
    // cout << "breakpoint A22. pars: " <<  range*10.0 << '\t' <<  0 << '\t' << range << endl;
    primaryTree->Project("primaryeSpectrPDGLogBinsHist", "ekin", TString::Format("pdg==%i || pdg==%i",pdg, -pdg));
    // cout << "breakpoint A3" << endl;
    primaryeSpectrPDGHist->SetLineColor(kRed);
    primaryeSpectrHist->GetXaxis()->SetRangeUser(lowRange, 1.5*highRange);
    primaryeSpectrPDGHist->GetXaxis()->SetRangeUser(lowRange, 1.5*highRange);
    if(verbosityLevel>0){
        TCanvas* primarySpectrumPDGCanvas = new TCanvas("primarySpectrumPDGCanvas", TString::Format("Energy spectrum of pdg==%i primaries", pdg));
        primarySpectrumPDGCanvas->SetGridx(); primarySpectrumPDGCanvas->SetGridy();
        primarySpectrumPDGCanvas->SetLogx(); primarySpectrumPDGCanvas->SetLogy();
        primaryeSpectrHist->DrawClone();
        primaryeSpectrPDGHist->DrawClone("same");
        cout << TString::Format("Total entries: %i, Chg. entries: %i\n", primaryeSpectrHist->GetEntries(), primaryeSpectrPDGHist->GetEntries());
    }
    // cout << "breakpoint B" << endl;
    std::vector<TH1D> retVar;
    retVar.push_back(*primaryeSpectrPDGHist);
    retVar.push_back(*primaryeSpectrHist);
    // cout << "breakpoint C" << endl;
    return retVar;
}
// Energy XY profile of the particles. Each point (x,y) is weightned with the corresponding energy
std::vector<TH2D> primaryEnergyProfile(){
    TH2D xyWeightnedHist = TH2D("xyWeightnedHist", "Energy profile of primaries in the transverse plane / BX; x [mm]; y [mm]; energy [GeV]", 200, -10., 10., 200, -10., 10.);
    xyWeightnedHist.GetXaxis()->SetRangeUser(-3,3);
    xyWeightnedHist.GetYaxis()->SetRangeUser(-3,3);

    // Initialize variables for the calculation
    double x0, y0, ekin;
    if(verbosityLevel>0) cout << "Setting branch addresses...";
    primaryTree->SetBranchAddress("x0", &x0);                       // mm
    if(verbosityLevel>0) cout << "x0,";
    primaryTree->SetBranchAddress("y0", &y0);   
    if(verbosityLevel>0) cout << "y0,";
    primaryTree->SetBranchAddress("ekin", &ekin);   
    if(verbosityLevel>0) cout << "ekin.\n";

    Long64_t primaryEntries = primaryTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< primaryEntries; i++){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);
        primaryTree->GetEntry(i);

        xyWeightnedHist.Fill(x0, y0, ekin);
    }

    if(verbosityLevel > 0){
        TCanvas* primaryEnergyProfileCanvas = new TCanvas("primaryEnergyProfileCanvas", "Energy XY profile of the particles");
        xyWeightnedHist.DrawClone("LEGO2Z");
    }
    std::vector<TH2D> retVar;
    retVar.push_back(xyWeightnedHist);
    return retVar;
}
// Event
// Energy spectrum of depositions for the upstream/downstream detector in a window Ea-Eb. Energy in keV
std::vector<TH1D> edepSpectrum(int bins, float emin, float emax){
    TH1D* totEdepdet0 = new TH1D("eSpectrumdet0", "Energy spectrum in the upstream detector;Deposited energy [keV];", bins, emin, emax);
    TH1D* totEdepdet1 = new TH1D("eSpectrumdet1", "Energy spectrum in the downstream detector;Deposited energy [keV];", bins, emin, emax);
    eventTree->Project("eSpectrumdet0", "edep", "det==0");
    eventTree->Project("eSpectrumdet1", "edep", "det==1");
    if(verbosityLevel>0){
        new TCanvas();
        totEdepdet0->DrawClone();
        new TCanvas();
        totEdepdet1->DrawClone();
    }

    std::vector<TH1D> retVar;
    retVar.push_back(*totEdepdet0);
    retVar.push_back(*totEdepdet1);
    delete totEdepdet0, totEdepdet1;
    return retVar;
}
// Energy spectrum of depositions for the upstream/downstream detector in a window Ea-Eb, partitioned by primary energies. Energy in keV 
std::vector<TH1D> edepSpectrumPartitionedByPrimaryEnergies(int bins=1000, double emin = 0.0 /*GeV*/, double emax= 6.0 /*GeV*/, double baseRatio = -1E6 /*keV*/){
    // Lambda function for formatting energy //input in GeV
    auto formatEnergy = [](double x){
        TString unit;
        if(x>=0.1){
            unit = TString::Format("%.2f GeV", x);
        }else if(0.001<=x && x<0.1){
            unit = TString::Format("%.2f MeV", x*1E3);
        }else if(0.000001<=x && x<0.001){
            unit = TString::Format("%.2f keV", x*1E6);
        }else if(x<0.000001){
            unit = TString::Format("%.2f eV", x*1E9);
        }else{
            unit = TString::Format("%f", x);
        }
        return unit;
    };

    // Variables
    double enLow=(emax+emin)/2;
    double enHigh=(emax+emin)/2;
    int slices;
    int eventIdPrimary; double ekinPrimary;
    int eventId, det; double edep;
    Long64_t primaryEntries, eventEntries;
    std::vector<double> slEnLow, slEnHigh;

    // Find min/max range of the primary energiesa and calculate the slices
    // Primary tree
    primaryTree->SetBranchAddress("event", &eventIdPrimary);
    primaryTree->SetBranchAddress("ekin", &ekinPrimary);
    primaryEntries = primaryTree->GetEntries();
    // Deposition tree
    eventTree->SetBranchAddress("event", &eventId);
    eventTree->SetBranchAddress("det", &det);
    eventTree->SetBranchAddress("edep", &edep);
    eventEntries = eventTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    if(verbosityLevel>0) cout << "eventTree has " << eventEntries << " entries. Looping over entries...\n";

    

    // Adjust selected range [emin, emax] according to the actual [enLow,enHigh] of the energies of primaries
    for(Long64_t i=0; i< primaryEntries; i++ ){
        primaryTree->GetEntry(i);
        if(emin < ekinPrimary && ekinPrimary < emax){
            if(enHigh < ekinPrimary) enHigh = ekinPrimary;
            if(enLow > ekinPrimary) enLow = ekinPrimary;
        }
    }
    if(verbosityLevel>0) cout << "Primary energy is between " << formatEnergy(enLow) << " < E < " << formatEnergy(enHigh) << endl;



    // Calculate how many slides are needed. If baseRatio > 0 then use a log slicing, else if baseRatio < 0 go linear with width = -baseRatio
    if(baseRatio < 0 ){
        // Linear slices. Width is fixed.
        if(verbosityLevel>0) cout << "Linear spaced slices. Slice width is: " << formatEnergy(-baseRatio*1E-6);
        slices = TMath::Nint(- (enHigh-enLow)*1E6 / baseRatio);
    }else if(baseRatio > 0 ){
        // Asking for Log spaced slices
        if(verbosityLevel>0) cout << "Log spaced slices. log("<<baseRatio<<", "<<(enHigh-enLow)*1E6<<")";
        slices = TMath::Nint(log(baseRatio, (enHigh-enLow)*1E6));
    }else{
        cout << "Invalid baseRatio parameter: " << baseRatio << endl;
    }
    if(verbosityLevel>0) cout << ". The number of slices is: " << slices << endl;
    
    //
    for(int i=0; i<slices; i++){
        double eLowi, eHighi; /*GeV*/
        if(baseRatio > 0 ){
            // Log slicing
            if(i==0){
            eLowi = enLow*1E6; /*GeV->keV*/
            eHighi = enLow*1E6 + pow(baseRatio,i+1);
            }else if(i == (slices-1)){
                eLowi = enLow*1E6 + pow(baseRatio,i);
                eHighi = enHigh*1E6;
            }else{
                eLowi = enLow*1E6 + pow(baseRatio,i);
                eHighi = enLow*1E6 + pow(baseRatio,i+1);
            }
        }else{
            // Linear slicing
            if(i==0){
            eLowi = enLow*1E6; /*GeV->keV*/
            eHighi = enLow*1E6 + (-baseRatio)*(i+1);
            }else if(i == (slices-1)){
                eLowi = enLow*1E6 + (-baseRatio)*(i);
                eHighi = enHigh*1E6;
            }else{
                eLowi = enLow*1E6 + (-baseRatio)*(i);
                eHighi = enLow*1E6 + (-baseRatio)*(i+1);
            }
        }
        
        slEnLow.push_back(eLowi*1E-6);
        slEnHigh.push_back(eHighi*1E-6);
    }
    

    // Define histogram variables
    if(verbosityLevel>0) cout << "Initializing histograms...";
    TH1D* edepStripSelHist[2][slices];
    for(int i=0; i<2; i++){
        for(int j=0; j<slices; j++){
            TString helper = formatEnergy(slEnLow[j]) + " < E < " + formatEnergy(slEnHigh[j]);
            // cout << helper << endl;
            edepStripSelHist[i][j] = new TH1D("edepStripSelHist"+whichDetShort(i)+TString::Format("%ith",j) , helper, bins, 0, 1000);
            edepStripSelHist[i][j]->GetXaxis()->SetTitle("energy [keV]");
            edepStripSelHist[i][j]->SetStats(0);
        }
    }
    if(verbosityLevel>0) cout << "OK!\n";

    // Loop over primaryTree, for each event understand in which slice they fill and then loop over entries of eventTree asmlong ad eventEDEP is equal to the one i am interested
    if(verbosityLevel>0) cout << "Loop over primaryTree...\n";
    int pslice0=0;
    int prevEventNb;
    int currentSlice=0;
    Long64_t lastJIndex = 1;
    int* counters = (int*)malloc((slices+1)*sizeof(int));
    int* counters2 = (int*)malloc((slices+1)*sizeof(int));
    for(int i=0; i< slices; i++){
        counters[i] = 0;
        counters2[i] = 0;
    }

    primaryTree->GetEntry(0);
    eventTree->GetEntry(0);
    resetProgress();
    for(Long64_t i=1; i< primaryEntries; i++ ){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);

        primaryTree->GetEntry(i);
        if( ekinPrimary < enLow || ekinPrimary > enHigh) continue;                               // Selects events in energy range

        currentSlice = (ekinPrimary-enLow)/(enHigh-enLow)*slices;
        if(currentSlice >= slices){
            // cout << currentSlice << endl;
            currentSlice--;
        }

        counters[currentSlice]++;                                                               // Counts numbers of events in the i-th energy window of the ith slice
        for(Long64_t j=lastJIndex; j<eventEntries; j++){
            eventTree->GetEntry(j);
            // // cout << eventIdPrimary << "\t" << eventId << "\t" << j << "\t" << det << "\t" << currentSlice << "\t" << ekinPrimary;

            if(eventId>eventIdPrimary){
                // cout << "\t\tbreak" << endl;
                lastJIndex = j;
                break;
            }else if(eventId == eventIdPrimary){
                // cout << "\t\tcase!" << endl;
                lastJIndex = j;
                edepStripSelHist[det][currentSlice]->Fill(edep);
                counters2[currentSlice]++;
            } 
        }
    }

    if(verbosityLevel>0){
        for(int slice=0; slice < slices; slice++){
            cout << "Numer of primaries with energy between " << TString::Format("[%f, %f] GeV are: ", emin+slice*(-baseRatio*1E-6), emin+(slice+1)*(-baseRatio*1E-6)) << counters[slice] << "\t" << counters2[slice] << endl;
        }
        for(int i=0; i<2; i++){
            TCanvas* edepSpectrumPartitionedCanvas = new TCanvas("edepSpectrumPartitioned"+whichDetShort(i)+"Canvas", "edepSpectrumPartitioned function");
            edepSpectrumPartitionedCanvas->SetGridx(); edepSpectrumPartitionedCanvas->SetGridy();
            double lSizeY;
            if(slices <= 8){
                lSizeY = 0.6;
            }else if(slices > 8 && slices < 40){
                lSizeY = 0.6 - 0.3*(slices-8)/32.0;
            }else{
                lSizeY = 0.;
            }
            TLegend* legend = new TLegend(1-0.5, lSizeY, 1-0.1, 1-0.1);
            legend->SetHeader("Initial beam particles energies accounted:","L"); // option "C" allows to center the header
            for(int j=0; j<slices; j++){
                edepStripSelHist[i][j]->SetLineColor(j+2);
                legend->AddEntry(edepStripSelHist[i][j], edepStripSelHist[i][j]->GetTitle() + TString::Format(" - %i entr.",(int)(edepStripSelHist[i][j]->GetEntries())),"l");
            }
            for(int j=0; j<slices; j++){
                edepStripSelHist[i][j]->SetStats(0);
                edepStripSelHist[i][j]->Draw("hist same");
            }
            legend->Draw();
            edepSpectrumPartitionedCanvas->Print("plots/edepSpectrumPartitioned"+whichDetShort(i)+".pdf");
        }
    }

    std::vector<TH1D> retVar;
    for(int i=0; i<2; i++){
        for(int j=0; j<slices; j++){
            retVar.push_back(*edepStripSelHist[i][j]);
        }
    }
    slEnLow.clear(); slEnHigh.clear();
    free(counters); free(counters2);
    return retVar;
}
// Impact of charged particles to the spectrum of deposited energy. Return a 2-array with the spectrum of energy deposited in upstream/downstream detector in the range [a,b] with #bins. Sel is a selector for the values \xi={"5.0", "7.0", "10.0"} 
std::vector<TH1D> edepSpectrumChg(int bins, double a, double b){
    // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
    TH1D* eSpectrChgHist[2];
    for (int i=0; i<2; i++){
        eSpectrChgHist[i] = new TH1D(TString::Format("edepSpectrChgDet%i",i), "Energy spectrum " + whichDetVerbose(i) + " detector", bins, a, b);
    }
    // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
    std::vector<int> pChgList;
    if(verbosityLevel>0) cout << "Allocated space is: " << pChgList.capacity() << endl;
    int eventSel, pdg;
    primaryTree->SetBranchAddress("event", &eventSel);
    primaryTree->SetBranchAddress("pdg", &pdg);
    Long64_t primaryEntries = primaryTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< primaryEntries; i++ ){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);
        
        primaryTree->GetEntry(i);
        if(pdg != 22) pChgList.push_back(i);
    }
    // cout << "Entry point. Bins, a, b " << bins << '\t' << a << '\t' << b << endl;
    if(verbosityLevel>0) cout << "Allocated space is: " << pChgList.capacity() << endl;
    Long64_t chgNb = pChgList.size();
    if(verbosityLevel>0) cout << "Charged primaries: " << chgNb << endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pChgList.begin(), pChgList.end());

    // Select energy deposition events with charged primaries
    int event, det; double edep;
    eventTree->SetBranchAddress("event", &event);
    eventTree->SetBranchAddress("det", &det);
    eventTree->SetBranchAddress("edep", &edep);
    Long64_t eventEntries = eventTree->GetEntries();
    if(verbosityLevel>0) cout << "eventTree has " << eventEntries << " entries. Looping over entries...\n";
    // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
    int p0 = 0;
    resetProgress();
    for(Long64_t i=0; i< eventEntries; i++ ){
        eventTree->GetEntry(i);
        float status = (float)(i+1) / eventEntries;
        printProgress(status);
        for(int i=p0; i<chgNb; i++){
            double val=pChgList[i];
            if(event < val){
                break;
            }else if(event == val){
                eSpectrChgHist[det]->Fill(edep);
                p0=i;
                break;
            }
        }
    }
    // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1;
    std::vector<TH1D> eSpectrTotHist = edepSpectrum(bins, a, b);
    verbosityLevel = verbosityLevel_bak;
    for(int i=0; i<2; i++){
        eSpectrChgHist[i]->GetXaxis()->SetTitle("Energy [keV]");
        eSpectrChgHist[i]->GetYaxis()->SetTitle("counts");

        eSpectrTotHist[i].GetXaxis()->SetTitle("Energy [keV]");
        eSpectrTotHist[i].GetYaxis()->SetTitle("counts");

        eSpectrTotHist[i].SetStats(0);
        eSpectrChgHist[i]->SetStats(0);
        eSpectrChgHist[i]->SetLineColor(kRed);
    }

    if(verbosityLevel>0){
        for(int i=0; i<2; i++){
            TCanvas* edepSpectrumChgCanvas = new TCanvas("edepSpectrumChg"+whichDetShort(i)+"Canvas", "Impact of charged particles to the spectrum of deposited energy");
            TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
            legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
            legend->AddEntry(&eSpectrTotHist[i],TString::Format("#gamma,e^{-},e^{+} - %i entries",(int)eSpectrTotHist[i].Integral()),"l");
            legend->AddEntry(eSpectrChgHist[i],TString::Format("e^{-},e^{+}    - %i entries",(int)eSpectrChgHist[i]->Integral()),"l");
            eSpectrTotHist[i].DrawClone("hist");
            eSpectrChgHist[i]->DrawClone("same hist");
            eSpectrTotHist[i].SetStats(1);
            eSpectrChgHist[i]->SetStats(1);
            legend->Draw();
        }
    }

    std::vector<TH1D> retVar;
    retVar.push_back(eSpectrTotHist[0]);
    retVar.push_back(eSpectrTotHist[1]);
    retVar.push_back(*eSpectrChgHist[0]);
    retVar.push_back(*eSpectrChgHist[1]);
    delete eSpectrChgHist[0], eSpectrChgHist[1];
    return retVar;
}
// Impact of charged particles to the spectrum of deposited energy. Default energy range 0-1000 and 1000 bins. Sel is a selector for the values \xi={"5", "7", "10"}
std::vector<TH1D> edepSpectrumChg(){
    return edepSpectrumChg(1000, 0, 1000);
}
// Histogram the energy spectrum of the hit stripNb
std::vector<TH1D> eSpectrumStrip(int stripNb = 100, float range = 5.){
    int binsNb = 100;
    TH1D* edepSpectrdet0Hist = new TH1D("edepSpectrdet0Hist", TString::Format("Energy spectrum for the strip %i for the upstream detector;energy [keV]; counts",stripNb), binsNb, 0., range);
    TH1D* edepSpectrdet1Hist = new TH1D("edepSpectrdet1Hist", TString::Format("Energy spectrum for the strip %i for the downstream detector;energy [keV]; counts",stripNb), binsNb, 0., range);
    for(int det=0; det<2; det++){
        detTree->Project(TString::Format("edepSpectrdet%iHist", det), "edep", TString::Format("det==%i && strip==%i", det, stripNb));
    }

    if(verbosityLevel>0){
        new TCanvas();
        edepSpectrdet0Hist->DrawClone("hist");
        new TCanvas();
        edepSpectrdet1Hist->DrawClone("hist");
    }

    std::vector<TH1D> retVar;
    retVar.push_back(*edepSpectrdet0Hist);
    retVar.push_back(*edepSpectrdet1Hist);
    delete edepSpectrdet0Hist, edepSpectrdet1Hist;
    return retVar;
    
}
// Return the total energy deposited in detNb in keV! Print out in GeV.
double totalEdep(int detNb=0){
    int det; double edepDet;
    eventTree->SetBranchAddress("det", &det);
    eventTree->SetBranchAddress("edep", &edepDet);
    Long64_t detEntries = eventTree->GetEntries();
    double totDetEn[2] = {0};
    for(Long64_t i=0; i< detEntries; i++ ){
            eventTree->GetEntry(i);
            if(det > 2) continue;
            totDetEn[det] += edepDet;
        }

    if(verbosityLevel > 0){
        cout << "File: " << ifilename << endl;
        cout << "Total energy deposited in the upstream detector is: " << totDetEn[0]/1E6 << " GeV" << endl;
        cout << "Total energy deposited in the downstream detector is: " << totDetEn[1]/1E6 << " GeV" << endl;
    }
    return totDetEn[detNb];     
}
// 2-array of energy deposited in the strip 100 for upstream/downstream detectors
double* eStrip_array(int stripNb = 100){
    // Total energy deposited in the strip
    int det, strip; double edep;
    detTree->SetBranchAddress("det", &det);               // 0,1
    detTree->SetBranchAddress("strip", &strip);           // 1-200
    detTree->SetBranchAddress("edep", &edep);             // keV
    static double totEnergyStrip[2] = {0};
    Long64_t detEntries = detTree->GetEntries();
    resetProgress();
    for(Long64_t i=0; i< detEntries; i++ ){
        detTree->GetEntry(i);
        float status = (float)(i+1)/(float)(detEntries);
        printProgress(status);

        if(strip != stripNb) continue;                                            // Select the correct detector
        totEnergyStrip[det] += edep;                                              // Total energy, useful for normalization (optional)
    }
    if(verbosityLevel>0){
        cout << "-----------------------------------------------------------------" << endl;
        cout << "File: " << ifilename << endl;
        cout << "-----------------------------------------------------------------" << endl;
        cout    << "Detector: upst. Energy deposited in strip #" << stripNb << " is: " << totEnergyStrip[0] << "keV" << endl;
        cout    << "Detector: down. Energy deposited in strip #" << stripNb << " is: " << totEnergyStrip[1] << "keV" << endl; 
    }
    return totEnergyStrip;
}
// Energy deposited in the strip 'stripNb' of detector 'detNb'
double eStrip(int detNb=0, int stripNb = 100){
    double result = eStrip_array(stripNb)[detNb];
    return result;
}
// Longitudinal energy deposition profile. Return 2-array with det0 and det1 plots
std::vector<TH1D> eLongitudinal_array(int nbBins = 100){
    TH1D* edepZHist[2];
    for(int i=0; i<2; i++){
        edepZHist[i] = new TH1D(TString::Format("edepZHistdet%i",i), "Longitudinal energy profile " + whichDetVerbose(i) + " detector;Z [mm]; energy [keV]", (int)nbBins, -0.050, 0.050);
        edepZHist[i]->GetXaxis()->SetRangeUser(-0.060, 0.060);
    }
    
    if(verbosityLevel>0) cout << "Filling the histogram..." << endl;
    int det; double z; double edep;
    debugTree->SetBranchAddress("det", &det);               //
    debugTree->SetBranchAddress("z", &z);                   // mm
    debugTree->SetBranchAddress("edep", &edep);             //

    Long64_t treeEntries = debugTree->GetEntries();
    resetProgress();
    for(Long64_t i=0; i< treeEntries; i++ ){
        float status = (float)(i+1)/(float)treeEntries;
        printProgress(status);

        debugTree->GetEntry(i);
        double relZ = z;
        if(det==1) relZ -= 20.0;
        edepZHist[det]->Fill(relZ, edep);
    }

    if(verbosityLevel>0){
        TCanvas* eLongitudinalCanvas = new TCanvas("eLongitudinalCanvas", "Longitudinal energy profile", 1200, 400);
        eLongitudinalCanvas->Divide(2);
        for(int i=0; i<2; i++){
            eLongitudinalCanvas->cd(i+1);
            edepZHist[i]->DrawClone("hist");
        }
    }

    std::vector<TH1D> retVar;
    retVar.push_back(*edepZHist[0]);
    retVar.push_back(*edepZHist[1]);
    delete edepZHist[0], edepZHist[1];
    return retVar;
}
// Longitudinal energy deposition profile. Return 2-array with det0 and det1 plots. Overloading 1
TH1D eLongitudinal(int detNb){
    return eLongitudinal_array(100)[detNb];
}
// Longitudinal energy deposition profile. Return 2-array with det0 and det1 plots. Overloading 2
TH1D eLongitudinal(int detNb, int nbBins){
    return eLongitudinal_array(nbBins)[detNb];
}
// Energy deposition in the transverse XY plane. Either X or Y profile by selecting the returned value [0] det0X [1] det0Y [2] det1X [3] det1Y
std::vector<TH1D> edepProject(double meshSize = 0.100){
    TH1D* result[4];
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = 0;
    std::vector<TH2D> temp = doseXY(meshSize);
    verbosityLevel = verbosityLevel_bak;
    for(int i=0; i<2; i++){
        TString helper = "Energy deposited in the " + whichDetVerbose(i) + " detector";
        result[2*i] = new TH1D(*temp[i].ProjectionX());
        result[2*i]->SetTitle(helper);
        result[2*i]->GetYaxis()->SetTitle("Deposited energy / BX [keV]");
        result[2*i+1] = new TH1D(*temp[i].ProjectionY());
        result[2*i+1]->SetTitle(helper);
        result[2*i+1]->GetYaxis()->SetTitle("Deposited energy / BX [keV]");
    }

    if(verbosityLevel>0){
        TCanvas* edepProjectCanvas = new TCanvas("edepProjectCanvas", "Energy deposition in the transverse XY plane", 1200, 600);
        edepProjectCanvas->Divide(2,2);
        for(int i=0; i<2; i++){
            edepProjectCanvas->cd(i+1);
            result[i]->DrawClone("hist");
            edepProjectCanvas->cd(i+3);
            result[i+2]->DrawClone("hist");
        }
    }

    std::vector<TH1D> retVar;
    for(int i=0; i<4; i++){
        retVar.push_back(*result[i]);
        delete result[i];
    }
    return retVar;
}
// Energy deposition profile in the strips. It uses the strip tree data. Return a 2-array for det0 and 1
std::vector<TH1D> edepProjectStripTree(){
    TH1D* result[2];
    for(int i=0; i<2; i++){
        TString helper = "Energy deposited in the " + whichDetVerbose(i) + " detector";
        result[i] = new TH1D("edepProjectStrip"+whichDetShort(i), helper, 200, 1, 200);
        result[i]->GetXaxis()->SetTitle("strip number [1-200]");
    }
    int det, strip; double edep;
    stripTree->SetBranchAddress("det", &det);
    stripTree->SetBranchAddress("strip", &strip);
    stripTree->SetBranchAddress("edep", &edep);
    Long64_t stripEntries = stripTree->GetEntries();
    if(verbosityLevel >0) cout << "stripTree has " << stripEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< stripEntries; i++ ){
        stripTree->GetEntry(i);
        float status = (float)(i+1) / stripEntries;
        printProgress(status);
        result[det]->Fill(strip, edep);
    }
    
    // Rescale to the BX
    int config = getConfig();
    switch(config){
        case 9:
        for(int i=0; i<2; i++){
            result[i]->Scale(1E-6);
            result[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 8:
        for(int i=0; i<2; i++){
            result[i]->Scale(1E-5);
            result[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 7:
        for(int i=0; i<2; i++){
            result[i]->Scale(1E-4);
            result[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 6:
        for(int i=0; i<2; i++){
            result[i]->Scale(1E-3);
            result[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        default:
        cout << "Config not found. getConfig() = " << config << endl;
    }

    if(verbosityLevel>0){
        TCanvas* edepProjectStripTreeCanvas = new TCanvas("edepProjectStripTreeCanvas", "Energy deposition profile in the strips", 1200, 400);
        edepProjectStripTreeCanvas->Divide(2);
        for(int i=0; i<2; i++){
            edepProjectStripTreeCanvas->cd(i+1);
            result[i]->DrawClone("hist");
        }
    }

    std::vector<TH1D> retVar;
    retVar.push_back(*result[0]);
    retVar.push_back(*result[1]);
    delete result[0], result[1];
    return retVar;
}
// Energy deposition profile in the strips given primaries in selected energy range [pE0, pE1] GeV. Return a 2-array for det0 and 1
std::vector<TH1D> edepProjectStripTree(double pE0, double pE1){
    // Get the list of tracks of primaries with energy in the given range
    std::vector<int> pSelList;
    int eventSel; double ekin;
    primaryTree->SetBranchAddress("event", &eventSel);
    primaryTree->SetBranchAddress("ekin", &ekin);
    Long64_t primaryEntries = primaryTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< primaryEntries; i++ ){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);
        
        primaryTree->GetEntry(i);
        if(pE0 < ekin && ekin < pE1) pSelList.push_back(i);
    }

    Long64_t selNb = pSelList.size();
    if(verbosityLevel>0) cout << "Numer of primaries with energy between " << TString::Format("[%f, %f] GeV: ", pE0, pE1) << selNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelList.begin(), pSelList.end());

    // Fill the strip energy depositions with selected events
    // Define variables
    TH1D* edepStripHist[2];
    TH1D* edepStripSelHist[2];
    for(int i=0; i<2; i++){
        TString helper = "Energy dep.s in strips for primaries with ";
        helper += TString::Format("E#in#left[%.2f, %.2f #right] GeV", pE0, pE1);
        helper += " for " + whichDetVerbose(i) + " detector";
        edepStripSelHist[i] = new TH1D("edepStripSelHist"+whichDetShort(i), helper, 200, 1, 200);
        edepStripHist[i] = new TH1D("edepStripHist"+whichDetShort(i), helper, 200, 1, 200);
        edepStripSelHist[i]->GetXaxis()->SetTitle("strip number [1-200]");
        edepStripHist[i]->GetXaxis()->SetTitle("strip number [1-200]");
    }
    // Fill histogram
    int event, det, strip; double edep;
    stripTree->SetBranchAddress("event", &event);
    stripTree->SetBranchAddress("det", &det);
    stripTree->SetBranchAddress("strip", &strip);
    stripTree->SetBranchAddress("edep", &edep);
    Long64_t stripEntries = stripTree->GetEntries();
    if(verbosityLevel >0) cout << "stripTree has " << stripEntries << " entries. Looping over entries...\n";
    
    int p0=0;
    resetProgress();
    for(Long64_t i=0; i< stripEntries; i++ ){
        stripTree->GetEntry(i);
        float status = (float)(i+1) / stripEntries;
        printProgress(status);

        edepStripHist[det]->Fill(strip, edep);

        for(int i=p0; i<selNb; i++){
            double val=pSelList[i];
            if(event < val){
                break;
            }else if(event == val){
                edepStripSelHist[det]->Fill(strip, edep);
                p0=i;
                break;
            }
        }
    }

    // Rescale to the BX
    int config = getConfig();
    switch(config){
        case 9:
        for(int i=0; i<2; i++){
            edepStripSelHist[i]->Scale(1E-6);
            edepStripSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepStripHist[i]->Scale(1E-6);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 8:
        for(int i=0; i<2; i++){
            edepStripSelHist[i]->Scale(1E-5);
            edepStripSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepStripHist[i]->Scale(1E-5);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 7:
        for(int i=0; i<2; i++){
            edepStripSelHist[i]->Scale(1E-4);
            edepStripSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepStripHist[i]->Scale(1E-4);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 6:
        for(int i=0; i<2; i++){
            edepStripSelHist[i]->Scale(1E-3);
            edepStripSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepStripHist[i]->Scale(1E-3);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        default:
        cout << "Config not found. getConfig() = " << config << endl;
    }

    std::vector<TH1D> retVar;
    for(int i=0; i<2; i++){
        retVar.push_back(*edepStripSelHist[i]);
        delete edepStripSelHist[i];
    }
    for(int i=0; i<2; i++){
        retVar.push_back(*edepStripHist[i]);
        delete edepStripHist[i];
    }
    return retVar;
}
// Energy deposition profile in the strips with comparison with deposition from charged primaries. It uses the strip tree data. Return a 2-array for det0 and 1
std::vector<TH1D> edepStripChg(){
    // Select events with the given energy range
    std::vector<int> pSelList;
    int eventSel, pdg;
    primaryTree->SetBranchAddress("event", &eventSel);
    primaryTree->SetBranchAddress("pdg", &pdg);
    Long64_t primaryEntries = primaryTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< primaryEntries; i++ ){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);
        
        primaryTree->GetEntry(i);
        if(pdg == -11 || pdg == 11) pSelList.push_back(i);
    }

    Long64_t selNb = pSelList.size();
    if(verbosityLevel>0) cout << "Charged primaries: " << selNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelList.begin(), pSelList.end());

    TH1D* edepChgHist[2];
    TH1D* edepStripHist[2];
    for(int i=0; i<2; i++){
        TString helper = "Energy deposited in strips of the " + whichDetVerbose(i) + " detector";
        edepChgHist[i] = new TH1D("edepChgHist"+whichDetShort(i), helper, 200, 1, 200);
        edepStripHist[i] = new TH1D("edepStripHist"+whichDetShort(i), helper, 200, 1, 200);
        edepChgHist[i]->GetXaxis()->SetTitle("strip number [1-200]");
        edepStripHist[i]->GetXaxis()->SetTitle("strip number [1-200]");
        edepChgHist[i]->SetLineColor(kRed);
    }
    int event, det, strip; double edep;
    stripTree->SetBranchAddress("event", &event);
    stripTree->SetBranchAddress("det", &det);
    stripTree->SetBranchAddress("strip", &strip);
    stripTree->SetBranchAddress("edep", &edep);
    Long64_t stripEntries = stripTree->GetEntries();
    if(verbosityLevel >0) cout << "stripTree has " << stripEntries << " entries. Looping over entries...\n";
    
    int p0=0;
    resetProgress();
    for(Long64_t i=0; i< stripEntries; i++ ){
        stripTree->GetEntry(i);
        float status = (float)(i+1) / stripEntries;
        printProgress(status);

        edepStripHist[det]->Fill(strip, edep);

        for(int i=p0; i<selNb; i++){
            double val=pSelList[i];
            if(event < val){
                break;
            }else if(event == val){
                edepChgHist[det]->Fill(strip, edep);
                p0=i;
                break;
            }
        }
    }

    // Rescale to the BX
    int config = getConfig();
    switch(config){
        case 9:
        for(int i=0; i<2; i++){
            edepStripHist[i]->Scale(1E-6);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepChgHist[i]->Scale(1E-6);
            edepChgHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 8:
        for(int i=0; i<2; i++){
            edepStripHist[i]->Scale(1E-5);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepChgHist[i]->Scale(1E-5);
            edepChgHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 7:
        for(int i=0; i<2; i++){
            edepStripHist[i]->Scale(1E-4);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepChgHist[i]->Scale(1E-4);
            edepChgHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        case 6:
        for(int i=0; i<2; i++){
            edepStripHist[i]->Scale(1E-3);
            edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            edepChgHist[i]->Scale(1E-3);
            edepChgHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        }
        break;

        default:
        cout << "Config not found. getConfig() = " << config << endl;
    }
    
    if(verbosityLevel>0){
        for(int i=0; i<2; i++){
            TCanvas* edepStripChgCanvas = new TCanvas("edepStripChg"+whichDetShort(i), "Energy deposited per strip in the " + whichDetVerbose(i) + " detector");
            TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
            edepStripHist[i]->DrawClone("hist");
            edepChgHist[i]->DrawClone("hist same");
            legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
            legend->AddEntry(edepStripHist[i], TString::Format("#gamma, others - %i entries",(int)edepStripHist[i]->GetEntries()),"l");
            legend->AddEntry(edepChgHist[i], TString::Format("e^{-},e^{+}    - %i entries",(int)edepChgHist[i]->GetEntries()),"l");
            legend->DrawClone();
            delete legend;
        }
    }

    std::vector<TH1D> retVar;
    retVar.push_back(*edepStripHist[0]);
    retVar.push_back(*edepStripHist[1]);
    retVar.push_back(*edepChgHist[0]);
    retVar.push_back(*edepChgHist[1]);
    delete edepStripHist[0], edepStripHist[1], edepChgHist[0], edepChgHist[1];
    return retVar;
}
// Energy depositions in the strips separated between the contributions with different energies
std::vector<TH1D> edepProjectPrimEnPartitioned(double enLow=0 /*GeV*/, double enHigh=10.0 /*GeV*/, double baseRatio = 10.0){
    // Lambda function for formatting
    auto formatEnergy = [](double x){
        TString unit;
        if(x>=0.1){
            unit = TString::Format("%.1f GeV", x);
        }else if(0.001<=x && x<0.1){
            unit = TString::Format("%.1f MeV", x*1E3);
        }else if(0.000001<=x && x<0.001){
            unit = TString::Format("%.1f keV", x*1E6);
        }else if(x<0.000001){
            unit = TString::Format("%.1f eV", x*1E9);
        }else{
            unit = TString::Format("%f", x);
        }
        return unit;
    };

    int verbosityLevel_bak = verbosityLevel; verbosityLevel = -1;
    const int slices = TMath::Nint(log(10.0, (enHigh-enLow)*1E6));
    double range = (enHigh-enLow)/slices;
    std::vector<TH1D> partitions;
    double maxCounts[2] = {0};
    if(verbosityLevel>0) cout << "The number of slice is: " << slices << endl;
    
    for(int i=0; i<slices; i++){
        double eLowi,eHighi;
        if(i==0){
            eLowi = enLow*1E6;
            eHighi = enLow*1E6 + pow(baseRatio,i+1);
        }else if(i == (slices-1)){
            eLowi = enLow*1E6 + pow(baseRatio,i);
            eHighi = enHigh*1E6;
        }else{
            eLowi = enLow*1E6 + pow(baseRatio,i);
            eHighi = enLow*1E6 + pow(baseRatio,i+1);
        }

        cout << "eLowi: " << formatEnergy(eLowi*1E-6) << " eHighi: " << formatEnergy(eHighi*1E-6) << endl;
        std::vector<TH1D> slice = edepProjectStripTree(eLowi*1E-6, eHighi*1E-6);
        for(int det=0; det<2; det++){
            TH1D tmpVar = slice[det];
            tmpVar.SetLineColor(2+i);
            tmpVar.SetTitle(formatEnergy(eLowi*1E-6) + " < E < " + formatEnergy(eHighi*1E-6) + TString::Format(" - %i entr.",(int)tmpVar.GetEntries()) );
            if(tmpVar.GetMaximum() > maxCounts[det]) maxCounts[det] = tmpVar.GetMaximum();
            partitions.push_back(tmpVar);
        }
    }

    std::vector<TH1D> fullHist = edepProjectStripTree();
    for(int det=0; det<2; det++){
        fullHist[det].SetLineColor(kBlue);
        fullHist[det].SetTitle(fullHist[det].GetTitle() + TString::Format(" - %i entr.",(int)fullHist[det].GetEntries()));
        partitions.push_back(fullHist[det]);
    }
    verbosityLevel = verbosityLevel_bak;

    // Debug
    if(verbosityLevel>0){
        for(int det=0; det<2; det++){
            TCanvas* edepSlicesCanvas = new TCanvas("edepSlices"+whichDetShort(det), "Report " + whichDetShort(det));
            edepSlicesCanvas->SetGridx(); edepSlicesCanvas->SetGridy();
            TLegend* legend = new TLegend(1.0-0.4, 0.7, 1.0-0.1, 0.9);
            fullHist[det].SetStats(0);
            fullHist[det].SetMinimum(1);
            fullHist[det].Draw("hist");
            legend->AddEntry(&fullHist[det], "0 < E < " + formatEnergy(fullHist[det].GetMaximum()), "l");
            for(int i= slices-1; i>=0; i--){
                // Make sure that all histograms are viewed correctly
                partitions[i*2+det].GetYaxis()->SetRangeUser(0, 1.05*maxCounts[det]);
                partitions[i*2+det].SetStats(0);
                partitions[i*2+det].SetMinimum(1);
                partitions[i*2+det].DrawClone("hist same");
                legend->AddEntry(&partitions[i*2+det], partitions[i*2+det].GetTitle(), "l");
            }
            legend->Draw();
        }
    }
    return partitions;
}
// Plot function edepProjectPrimEnPartitioned
void plot_edepProjectPrimEnPartitioned(double enLow=0 /*GeV*/, double enHigh=10.0 /*GeV*/, double baseRatio = 10.0, TString path = "plots/", TString format = ".pdf"){
    // Lambda function for formatting
    auto formatEnergy = [](double x){
        TString unit;
        if(x>=0.1){
            unit = TString::Format("%.1f GeV", x);
        }else if(0.001<=x && x<0.1){
            unit = TString::Format("%.1f MeV", x*1E3);
        }else if(0.000001<=x && x<0.001){
            unit = TString::Format("%.1f keV", x*1E6);
        }else if(x<0.000001){
            unit = TString::Format("%.1f eV", x*1E9);
        }else{
            unit = TString::Format("%f", x);
        }
        return unit;
    };

    int verbosityLevel_bak = verbosityLevel; verbosityLevel = -1;
    const int slices = TMath::Nint(log(10.0, (enHigh-enLow)*1E6));
    double range = (enHigh-enLow)/slices;
    TH1D partitions[slices][2];
    double maxCounts[2] = {0};
    if(verbosityLevel>0) cout << "The number of slice is: " << slices << endl;
    std::vector<TH1D> fullHist = edepProjectStripTree();
    for(int i=0; i<slices; i++){
        double eLowi,eHighi;
        if(i==0){
            eLowi = enLow*1E6;
            eHighi = enLow*1E6 + pow(baseRatio,i+1);
        }else if(i == (slices-1)){
            eLowi = enLow*1E6 + pow(baseRatio,i);
            eHighi = enHigh*1E6;
        }else{
            eLowi = enLow*1E6 + pow(baseRatio,i);
            eHighi = enLow*1E6 + pow(baseRatio,i+1);
        }

        cout << "eLowi: " << formatEnergy(eLowi*1E-6) << " eHighi: " << formatEnergy(eHighi*1E-6) << endl;
        std::vector<TH1D> temp = edepProjectStripTree(eLowi*1E-6, eHighi*1E-6);
        for(int j=0; j<2; j++){
            partitions[i][j] = temp[j];
        }
        
        for(int det=0; det<2; det++){
            fullHist[det].SetLineColor(kBlue);
            partitions[i][det].SetLineColor(2+i);
            partitions[i][det].SetTitle(formatEnergy(eLowi*1E-6) + " < E < " + formatEnergy(eHighi*1E-6));
            if(partitions[i][det].GetMaximum() > maxCounts[det]) maxCounts[det] = partitions[i][det].GetMaximum();
        }
        // new TCanvas();
        // partitions[i][0].DrawClone("hist same");
    }

    for(int det=0; det<2; det++){
        TCanvas painter("edepSlices"+whichDetShort(det), "Report " + whichDetShort(det));
        TLegend* legend = new TLegend(1.0-0.4, 0.7, 1.0-0.1, 0.9);
        fullHist[det].SetStats(0);
        fullHist[det].SetMinimum(1);
        fullHist[det].Draw("hist");
        legend->AddEntry(&fullHist[det], "0 < E < " + formatEnergy(fullHist[det].GetMaximum()), "l");
        for(int i= slices-1; i>=0; i--){
            // Make sure that all histograms are viewed correctly
            partitions[i][det].GetYaxis()->SetRangeUser(0, 1.05*maxCounts[det]);
            partitions[i][det].SetStats(0);
            partitions[i][det].SetMinimum(1);
            partitions[i][det].Draw("hist same");
            legend->AddEntry(&partitions[i][det], partitions[i][det].GetTitle(), "l");
        }
        legend->Draw();
        painter.SetGridx(); painter.SetGridy();
        TString fname = painter.GetName();
        painter.Print(path+fname+format);
        painter.SetLogy();
        painter.Print(path+fname+"Log"+format);
    }

    verbosityLevel = verbosityLevel_bak;
}


/*
* ******************************************************************************************************************************************
* Dose
* ******************************************************************************************************************************************
*/
// Energy & dose 2D map using a mesh of size LxLx100 um3. Return a 4-array of TH2D with [0]emapdet0, [1]emapdet1, [2]dmapdet0, [3]dmapdet1
std::vector<TH2D> doseXY(double meshSize = 0.100 /*mm*/){
    const int nbXpts = 20./meshSize;
    const int nbYpts = nbXpts;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector [keV];X [mm]; Y [mm]", nbXpts, -10, 10, nbYpts, -10, 10);
        // edepXYHist[i].GetXaxis()->SetRangeUser(-1.5, 1.5);
        // edepXYHist[i].GetYaxis()->SetRangeUser(-1.5, 1.5);
        // Style
        // edepXYHist[i]->SetStats(0);
        edepXYHist[i]->SetContour(500);
    }

    if(verbosityLevel>0) cout << "Filling the histogram..." << endl;
    int det, pdg; double x, y; double edep;
    debugTree->SetBranchAddress("det", &det);               // 0,1
    debugTree->SetBranchAddress("x", &x);                   // mm
    debugTree->SetBranchAddress("y", &y);                   // mm
    debugTree->SetBranchAddress("edep", &edep);             // keV
    debugTree->SetBranchAddress("pdg", &pdg);               // -11, 11, 22, -

    Long64_t treeEntries = debugTree->GetEntries();
    resetProgress();
    for(Long64_t i=0; i< treeEntries; i++ ){
        debugTree->GetEntry(i);
        float status = (float)(i+1)/(float)treeEntries;
        printProgress(status);
        if(det >= 2) continue;
        edepXYHist[det]->Fill(x, y, edep);
        totEnergy[det] += edep;
    }

    double conversionFactor = 0.1 * 1.60 / 3.970;
    double vol = (meshSize*1E3)*(meshSize*1E3)*(100);
    TH2D* doseXYHist[2];
    for(int i=0; i<2; i++){
        doseXYHist[i] = new TH2D(*edepXYHist[i]);
        doseXYHist[i]->SetName(TString::Format("doseXYHist%i",i));
        doseXYHist[i]->SetTitle("Dose XY map in the " + whichDetVerbose(i) + " detector [Gy] ;X [mm]; Y [mm]");
        doseXYHist[i]->Scale(conversionFactor/vol);
        // doseXYHist->GetXaxis()->SetRangeUser(-1.5, 1.5);
        // doseXYHist->GetYaxis()->SetRangeUser(-1.5, 1.5);
    }

    // Styling
    // gStyle->SetPalette(kRainBow);
    if(verbosityLevel>0){
        TCanvas* doseXYCanvas = new TCanvas("doseXYCanvas", "Energy & dose 2D map using a mesh of size LxLx100 um3", 1200, 800);
        doseXYCanvas->Divide(2,2);
        doseXYCanvas->SetFrameBorderSize(2);
        for(int i=0; i<2; i++){
            doseXYCanvas->cd(1+i);
            edepXYHist[i]->DrawClone("colz");
            doseXYCanvas->cd(3+i);
            doseXYHist[i]->DrawClone("colz");
        }
    }

    // Verbose messages
    if(verbosityLevel>1){
        cout << "Mesh size [um]: "  << meshSize*1000.           << endl;
        cout << "Mesh blocks: "     << nbXpts<<'x'<<nbYpts      << endl;
        cout    << "\n\n";
    }

    // Statistics & output
    double radius = gausBSigma[0]*sqrt(2.*log(20.));
    double totVol = TMath::Pi() * radius*radius * detThickness;
    double peakEnergy[2], peakDose[2], peakPoint[2][2];
    for(int i=0; i<2; i++){
        peakEnergy[i] = edepXYHist[i]->GetMaximum();
        peakDose[i] = doseXYHist[i]->GetMaximum();
        peakPoint[i][0] = doseXYHist[i]->GetXaxis()->GetBinCenter(doseXYHist[i]->GetMaximumBin());
        peakPoint[i][1] = doseXYHist[i]->GetYaxis()->GetBinCenter(doseXYHist[i]->GetMaximumBin());

        if(verbosityLevel>0){
            // Report
            cout    << "Detector: "                                 << whichDetVerbose(i) << endl;
            if(ifilename.Contains("_7_")){
                if(verbosityLevel>1){
                    cout    << "Total energy deposit is:\t"             << totEnergy[i]/1E6 << " GeV" << endl;
                    cout    << "Total absorbed dose is:\t\t"            << totEnergy[i]/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
                    cout    << "Peak dose is:\t\t\t"                    << peakDose[i] << " Gy (" << peakEnergy[i] << " keV in " << vol << " um3)" <<  endl;
                    cout    << "----------------BX--------------------" << endl;
                }
                cout    << "Total energy deposit/BX is:\t"          << totEnergy[i]/1E6 * 1E2 << " GeV" << endl;
                cout    << "Total absorbed dose/BX is:\t"           << totEnergy[i]/totVol * conversionFactor * 1E2 << " Gy (" << totVol << " um3)" << endl;
                cout    << "Peak dose/BX is:\t\t"                   << peakDose[i] * 1E2 << " Gy (" << peakEnergy[i] * 1E2 << " keV in " << vol << " um3)" <<  endl;
            }else if(ifilename.Contains("_9_")){    
                cout    << "Total energy deposit/BX is:\t"          << totEnergy[i]/1E6 << " GeV" << endl;
                cout    << "Total absorbed dose/BX is:\t"           << totEnergy[i]/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
                cout    << "Peak dose/BX is:\t\t"                   << peakDose[i] << " Gy (" << peakEnergy[i] << " keV in " << vol << " um3)" <<  endl;
            }
            if(i!=1)    cout    << "--------------------------------------" << endl;
        }
    }


    std::vector<TH2D> retVar;
    retVar.push_back(*edepXYHist[0]);
    retVar.push_back(*edepXYHist[1]);
    retVar.push_back(*doseXYHist[0]);
    retVar.push_back(*doseXYHist[1]);
    for(int i=0; i<2; i++){
        delete edepXYHist[i], doseXYHist[i];
    }
    return retVar;
}
// Energy & dose 2D map using a mesh of size LxLx100 um3. Contributions from only photons to the dose. Return a 4-array of TH2D with [0]emapdet0, [1]emapdet1, [2]dmapdet0, [3]dmapdet1
std::vector<TH2D> doseXY(bool neutral, double meshSize = 0.100 /*mm*/){
    // Select events with the given energy range
    std::vector<int> pSelList;
    int eventSel, pdgSel;
    primaryTree->SetBranchAddress("event", &eventSel);
    primaryTree->SetBranchAddress("pdg", &pdgSel);
    Long64_t primaryEntries = primaryTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< primaryEntries; i++ ){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);
        
        primaryTree->GetEntry(i);
        if(pdgSel != 2212) pSelList.push_back(i);
    }

    Long64_t selNb = pSelList.size();
    if(verbosityLevel>0) cout << "Photons in the primaries: " << selNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelList.begin(), pSelList.end());
    



    const int nbXpts = 20./meshSize;
    const int nbYpts = nbXpts;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector [keV];X [mm]; Y [mm]", nbXpts, -10, 10, nbYpts, -10, 10);
        // edepXYHist[i].GetXaxis()->SetRangeUser(-1.5, 1.5);
        // edepXYHist[i].GetYaxis()->SetRangeUser(-1.5, 1.5);
        // Style
        // edepXYHist[i]->SetStats(0);
        edepXYHist[i]->SetContour(500);
    }

    if(verbosityLevel>0) cout << "Filling the histogram..." << endl;
    int event, det, pdg; double x, y; double edep;
    debugTree->SetBranchAddress("event", &event);
    debugTree->SetBranchAddress("det", &det);               // 0,1
    debugTree->SetBranchAddress("x", &x);                   // mm
    debugTree->SetBranchAddress("y", &y);                   // mm
    debugTree->SetBranchAddress("edep", &edep);             // keV
    debugTree->SetBranchAddress("pdg", &pdg);               // -11, 11, 22, -
    Long64_t treeEntries = debugTree->GetEntries();

    int p0=0;
    resetProgress();
    for(Long64_t i=0; i< treeEntries; i++ ){
        debugTree->GetEntry(i);
        float status = (float)(i+1)/(float)treeEntries;
        printProgress(status);
        if(det >= 2) continue;
        totEnergy[det] += edep;

        for(int j=p0; j<selNb; j++){
            int val=pSelList[j];
            if(event < val){
                break;
            }else if(event == val){
                edepXYHist[det]->Fill(x, y, edep);
                p0=j;
                break;
            }
        }
    }

    double conversionFactor = 0.1 * 1.60 / 3.970;
    double vol = (meshSize*1E3)*(meshSize*1E3)*(100);
    TH2D* doseXYHist[2];
    for(int i=0; i<2; i++){
        doseXYHist[i] = new TH2D(*edepXYHist[i]);
        doseXYHist[i]->SetName(TString::Format("doseXYHist%i",i));
        doseXYHist[i]->SetTitle("Dose XY map in the " + whichDetVerbose(i) + " detector [Gy] ;X [mm]; Y [mm]");
        doseXYHist[i]->Scale(conversionFactor/vol);
        // doseXYHist->GetXaxis()->SetRangeUser(-1.5, 1.5);
        // doseXYHist->GetYaxis()->SetRangeUser(-1.5, 1.5);
    }

    // Styling
    // gStyle->SetPalette(kRainBow);
    if(verbosityLevel>0){
        TCanvas* doseXYCanvas = new TCanvas("doseXYCanvas", "Energy & dose 2D map using a mesh of size LxLx100 um3", 1200, 800);
        doseXYCanvas->Divide(2,2);
        doseXYCanvas->SetFrameBorderSize(2);
        for(int i=0; i<2; i++){
            doseXYCanvas->cd(1+i);
            edepXYHist[i]->DrawClone("colz");
            doseXYCanvas->cd(3+i);
            doseXYHist[i]->DrawClone("colz");
        }
    }

    // Verbose messages
    if(verbosityLevel>1){
        cout << "Mesh size [um]: "  << meshSize*1000.           << endl;
        cout << "Mesh blocks: "     << nbXpts<<'x'<<nbYpts      << endl;
        cout    << "\n\n";
    }

    // Statistics & output
    double radius = gausBSigma[0]*sqrt(2.*log(20.));
    double totVol = TMath::Pi() * radius*radius * detThickness;
    double peakEnergy[2], peakDose[2], peakPoint[2][2];
    for(int i=0; i<2; i++){
        peakEnergy[i] = edepXYHist[i]->GetMaximum();
        peakDose[i] = doseXYHist[i]->GetMaximum();
        peakPoint[i][0] = doseXYHist[i]->GetXaxis()->GetBinCenter(doseXYHist[i]->GetMaximumBin());
        peakPoint[i][1] = doseXYHist[i]->GetYaxis()->GetBinCenter(doseXYHist[i]->GetMaximumBin());

        if(verbosityLevel>0){
            // Report
            cout    << "Detector: "                                 << whichDetVerbose(i) << endl;
            if(ifilename.Contains("_7_")){
                if(verbosityLevel>1){
                    cout    << "Total energy deposit is:\t"             << totEnergy[i]/1E6 << " GeV" << endl;
                    cout    << "Total absorbed dose is:\t\t"            << totEnergy[i]/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
                    cout    << "Peak dose is:\t\t\t"                    << peakDose[i] << " Gy (" << peakEnergy[i] << " keV in " << vol << " um3)" <<  endl;
                    cout    << "----------------BX--------------------" << endl;
                }
                cout    << "Total energy deposit/BX is:\t"          << totEnergy[i]/1E6 * 1E2 << " GeV" << endl;
                cout    << "Total absorbed dose/BX is:\t"           << totEnergy[i]/totVol * conversionFactor * 1E2 << " Gy (" << totVol << " um3)" << endl;
                cout    << "Peak dose/BX is:\t\t"                   << peakDose[i] * 1E2 << " Gy (" << peakEnergy[i] * 1E2 << " keV in " << vol << " um3)" <<  endl;
            }else if(ifilename.Contains("_9_")){    
                cout    << "Total energy deposit/BX is:\t"          << totEnergy[i]/1E6 << " GeV" << endl;
                cout    << "Total absorbed dose/BX is:\t"           << totEnergy[i]/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
                cout    << "Peak dose/BX is:\t\t"                   << peakDose[i] << " Gy (" << peakEnergy[i] << " keV in " << vol << " um3)" <<  endl;
            }
            if(i!=1)    cout    << "--------------------------------------" << endl;
        }
    }


    std::vector<TH2D> retVar;
    retVar.push_back(*edepXYHist[0]);
    retVar.push_back(*edepXYHist[1]);
    retVar.push_back(*doseXYHist[0]);
    retVar.push_back(*doseXYHist[1]);
    for(int i=0; i<2; i++){
        delete edepXYHist[i], doseXYHist[i];
    }
    return retVar;
}
//NEW
// Energy deposition in the transverse XY plane. Either X or Y profile by selecting the returned value [0] det0X [1] det0Y [2] det1X [3] det1Y
std::vector<TH1D> edepProject(bool filter, double meshSize = 0.100){
    // Select events with the given filter (for now hard-coded below). Now is selects non-gamma particles
    // Upstream detector.
    // Selection happens using the event ID since primary particles are spawned at the surface of the upstream sapphire.
    std::vector<int> pSelListU;
    int eventSelU, pdgSelU;
    primaryTree->SetBranchAddress("event", &eventSelU);
    primaryTree->SetBranchAddress("pdg", &pdgSelU);
    Long64_t primaryEntries = primaryTree->GetEntries();
    if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< primaryEntries; i++ ){
        float status = (float)(i+1) / primaryEntries;
        printProgress(status);
        
        primaryTree->GetEntry(i);
        if(pdgSelU != 22) pSelListU.push_back(i);
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the upstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListU.begin(), pSelListU.end());
    
    // Downstream detector.
    // Selection happens using the track ID and the charge scoring plane tree
    class pSelListDCl {
        public:
        pSelListDCl(int event, int track){ eventID = event; trackID= track;};
        ~pSelListDCl(){};

        int getEvent(){ return eventID;}
        int getTrack(){ return trackID;}

        private:
        int eventID;
        int trackID;
    };
    std::vector<pSelListDCl> pSelListD; // Track IDs of the particles entering the downstream detector, and satifying the filter condition
    int eventSelD, trackSelD, pdgSelD;
    planeTree->SetBranchAddress("event", &eventSelD);
    planeTree->SetBranchAddress("track", &trackSelD);
    planeTree->SetBranchAddress("pdg", &pdgSelD);
    Long64_t planeEntries = planeTree->GetEntries();
    if(verbosityLevel>0) cout << "planeTree has " << planeEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< planeEntries; i++ ){
        float status = (float)(i+1) / planeEntries;
        printProgress(status);
        
        planeTree->GetEntry(i);
        if(pdgSelD != 22){
            pSelListDCl tmp(eventSelD, trackSelD);
            pSelListD.push_back(tmp);
        }
    }
    Long64_t selDNb = pSelListD.size();
    if(verbosityLevel>0) cout << "Non-photons entering the downstream detector: " << selDNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    //std::sort(pSelListD.begin(), pSelListD.end());    // They are sorted by contruction

    // Dose code with some minor modifications
    const int binsNb = 20./meshSize;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    TH2D* edepXYSelHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector [keV];X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        edepXYSelHist[i] = new TH2D(TString::Format("edepXYSelHist%i",i), "non-photons primaries " + whichDetVerbose(i) + " detector [keV];X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        // Style
        // edepXYHist[i]->SetStats(0);
        edepXYHist[i]->SetContour(500);
        edepXYSelHist[i]->SetContour(500);
        edepXYSelHist[i]->SetLineColor(kRed);
    }

    if(verbosityLevel>0) cout << "Filling the histogram..." << endl;
    int event, track, det, pdg; double x, y; double edep;
    debugTree->SetBranchAddress("event", &event);           // Event ID
    debugTree->SetBranchAddress("track", &track);           // Track ID
    debugTree->SetBranchAddress("det", &det);               // 0,1
    debugTree->SetBranchAddress("x", &x);                   // mm
    debugTree->SetBranchAddress("y", &y);                   // mm
    debugTree->SetBranchAddress("edep", &edep);             // keV
    debugTree->SetBranchAddress("pdg", &pdg);               // -11, 11, 22, -
    Long64_t debugEntries = debugTree->GetEntries();

    int pU0=0; int pD0=0;
    resetProgress();
    for(Long64_t i=0; i< debugEntries; i++ ){
        debugTree->GetEntry(i);
        float status = (float)(i+1)/(float)debugEntries;
        printProgress(status);
        if(det >= 2) continue;  
        totEnergy[det] += edep;
        // Filling energy depositions from all the primaries.
        edepXYHist[det]->Fill(x, y, edep);

        // Filling energy depositions from selected 'filter' primaries
        if(det==0){
            for(int j=pU0; j<selUNb; j++){
                int val=pSelListU[j];
                if(event < val){
                    break;
                }else if(event == val){
                    // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\tchg\n";
                    edepXYSelHist[det]->Fill(x, y, edep);
                    pU0=j;
                    break;
                }
            }
        }else{
            for(int j_events=pD0; j_events<selDNb; j_events++){
                auto pSelClass_ith = pSelListD[j_events];
                int evtOf_ith_class = pSelClass_ith.getEvent();
                int trkOf_ith_class = pSelClass_ith.getTrack();
                if(event < evtOf_ith_class){
                    break;
                }else if(event == evtOf_ith_class){
                    if(track == trkOf_ith_class){
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\tchg\n";
                        edepXYSelHist[det]->Fill(x, y, edep);
                    }
                    pD0=j_events;
                    break;
                }
            }
        }
    }


    // Verbose messages
    if(verbosityLevel>1){
        cout << "Mesh size [um]: "  << meshSize*1000.           << endl;
        cout << "Mesh blocks: "     << binsNb<<'x'<<binsNb      << endl;
        cout    << "\n\n";
    }
    
    std::vector<TH1D> edepProj;
    for(int i=0; i<2; i++){
        edepProj.push_back(TH1D(*edepXYHist[i]->ProjectionX()));
        edepProj.push_back(TH1D(*edepXYHist[i]->ProjectionY()));
    }
    for(int i=0; i<2; i++){
        edepProj.push_back(TH1D(*edepXYSelHist[i]->ProjectionX()));
        edepProj.push_back(TH1D(*edepXYSelHist[i]->ProjectionY()));
    }
    for(int i=0; i<8; i++){
        edepProj[i].GetYaxis()->SetTitle("Deposited energy / BX [keV]");
    }

    if(verbosityLevel>0){
        TCanvas* edepProjectCanvas = new TCanvas("edepProjectCanvas", "Energy deposition in the transverse XY plane", 1200, 600);
        edepProjectCanvas->Divide(2,2);
        for(int j=0; j<2; j++){
            for(int i=2*j; i<2*j+2; i++){
                edepProjectCanvas->cd(i+1);
                edepProjectCanvas->SetGridx(); edepProjectCanvas->SetGridy();
                TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
                legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
                legend->AddEntry(&edepProj[i],TString::Format("#gamma, other - %i entries",(int)edepProj[i].GetEntries()),"l");
                legend->AddEntry(&edepProj[i+4],TString::Format("other    - %i entries",(int)edepProj[i+4].GetEntries()),"l");
                edepProj[i].DrawClone("hist");
                edepProj[i+4].DrawClone("hist same");
                legend->Draw();
            }
        }
    }
    
    return edepProj;
}


/*
* ******************************************************************************************************************************************
* Detector performance for TDR - START
* ******************************************************************************************************************************************
*/
// Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Threshold is in keV. Return a 2-array with fit for upstream and downstream det.s
std::vector<TH1I> fitHits(double threshold, int* fitRange){
    static double rchisq[2];
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1; 
    std::vector<TH1D> edepStrip = edepProject();
    verbosityLevel = verbosityLevel_bak;
    TH1I* hitHist[2];
    for(int i=0; i<2; i++){
        hitHist[i] = new TH1I(TString::Format("hitXdet%i", i), "Hits in the "+whichDetVerbose(i)+" detector;Strip number [1-200];counts", 200, 1, 200);
    }

    // Initialize variables for the calculation
    int det, strip; double edep;
    if(verbosityLevel>0) cout << "Setting branch addresses...";
    stripTree->SetBranchAddress("det", &det);                  // 0,1
    if(verbosityLevel>0) cout << "det,";
    stripTree->SetBranchAddress("strip", &strip);              // 0-200
    if(verbosityLevel>0) cout << "strip,";
    stripTree->SetBranchAddress("edep", &edep);
    if(verbosityLevel>0) cout << "edep.\n";
    Long64_t stripEntries = stripTree->GetEntries();
    if(verbosityLevel>0) cout << "stripTree has " << stripEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< stripEntries; i++ ){
        float status = (float)(i+1) / stripEntries;
        printProgress(status);

        stripTree->GetEntry(i);
        if(edep < threshold) continue;                                         // Hit condition. If edep is less than threshold keV then quit the loop istance
        hitHist[det]->Fill(strip, (int)(edep/threshold));                      // Add the number of hits in the corresponding x/y det hist.
                                                                               // For example: edep=5.4 keV -> 5 hits
    }
    // I've noticed that error are not properly calculated, despite the option is
    //      root [13] result->GetBinErrorOption()
    //      (TH1::EBinErrorOpt) (TH1::kNormal) : (int) 0)
    // correct. So I wrote the following code to fix the issue
    if(verbosityLevel>0) cout << "Fixing errors...\n";
    resetProgress();
    for(int n=0; n < 2 ; n++ ){
        Long64_t binNb = hitHist[n]->GetNbinsX();
        for(Long64_t bin=0; bin < binNb; bin++){
            float status = (float)(bin*(n+1)+1)/(float)(2 * binNb);
            printProgress(status);
            double binContent = hitHist[n]->GetBinContent(bin);
            hitHist[n]->SetBinError(bin, sqrt(binContent));
        }
    }

    // Rescale to the BX
    int config = getConfig();
    switch(config){
        case 9:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-6);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i]->GetYaxis()->SetTitle("hits");
        }
        break;

        case 8:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-5);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i]->Scale(1E1);
            hitHist[i]->GetYaxis()->SetTitle("hits");
        }
        break;

        case 7:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-4);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i]->Scale(1E2);
            hitHist[i]->GetYaxis()->SetTitle("hits");
        }
        break;

        case 6:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-3);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i]->Scale(1E2);
            hitHist[i]->GetYaxis()->SetTitle("hits");
        }
        break;

        default:
        cout << "Config not found. getConfig() = " << config << endl;
    }

    // Fitting
    gStyle->SetOptFit(1);
    TString verbQuiet = "";
    if(verbosityLevel < 0) verbQuiet = "0Q";
    for(Int_t i=0; i<2; i++){
        TF1 fit = TF1("fit", "gaus", fitRange[0], fitRange[1]);
        fit.SetParameter(0, hitHist[i]->GetMaximum());
        fit.SetParameter(1, hitHist[i]->GetMean());
        fit.SetParameter(2, hitHist[i]->GetRMS());
        hitHist[i]->Fit("fit", "R"+verbQuiet);
        rchisq[i] = fit.GetChisquare();
    }

    if(verbosityLevel > 0){
        // Output plots
        TCanvas* fitHitsCanvas = new TCanvas("fitHitsCanvas", "Plots", 1000, 600);
        fitHitsCanvas->DivideSquare(4);
        for(int i=0; i<2; i++){
            fitHitsCanvas->cd(3+i);
            hitHist[i]->DrawClone("e");
        }
    }

    std::vector<TH1I> retVar;
    for(int i=0; i<2; i++){
        retVar.push_back(*hitHist[i]);
        delete hitHist[i];
    }
    return retVar;
}
// Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Overloading
std::vector<TH1I> fitHits(){
    int defRange[2] = {90, 110};
    return fitHits(1.0, defRange);
}
// Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Threshold is in keV, fitRange is a 2-array with stripL, stripH where the fit is ranged.
double* fitHitsChi2(double threshold, int* fitRange){
    static double rchisq[2];
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1; 
    std::vector<TH1D> edepStrip = edepProject();
    verbosityLevel = verbosityLevel_bak;
    TH1I hitHist[2];
    hitHist[0] = TH1I("hitXdet0", "Hits in the upstream detector;Strip number [1-200];counts", 200, 1, 200);
    hitHist[1] = TH1I("hitYdet1", "Hits in the downstream detector;Strip number [1-200];counts", 200, 1, 200);

    // Initialize variables for the calculation
    int det, strip; double edep;
    if(verbosityLevel>0) cout << "Setting branch addresses...";
    stripTree->SetBranchAddress("det", &det);                  // 0,1
    if(verbosityLevel>0) cout << "det,";
    stripTree->SetBranchAddress("strip", &strip);              // 0-200
    if(verbosityLevel>0) cout << "strip,";
    stripTree->SetBranchAddress("edep", &edep);
    if(verbosityLevel>0) cout << "edep.\n";
    Long64_t stripEntries = stripTree->GetEntries();
    if(verbosityLevel>0) cout << "stripTree has " << stripEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< stripEntries; i++ ){
        float status = (float)(i+1) / stripEntries;
        printProgress(status);

        stripTree->GetEntry(i);
        if(edep < threshold) continue;                                         // Hit condition. If edep is less than threshold keV then quit the loop istance
        hitHist[det].Fill(strip, (int)(edep/threshold));                      // Add the number of hits in the corresponding x/y det hist.
                                                                               // For example: edep=5.4 keV -> 5 hits
    }
    // I've noticed that error are not properly calculated, despite the option is
    //      root [13] result->GetBinErrorOption()
    //      (TH1::EBinErrorOpt) (TH1::kNormal) : (int) 0)
    // correct. So I wrote the following code to fix the issue
    if(verbosityLevel>0) cout << "Fixing errors...\n";
    resetProgress();
    for(int n=0; n < 2 ; n++ ){
        Long64_t binNb = hitHist[n].GetNbinsX();
        for(Long64_t bin=0; bin < binNb; bin++){
            float status = (float)(bin*(n+1)+1)/(float)(2 * binNb);
            printProgress(status);
            double binContent = hitHist[n].GetBinContent(bin);
            hitHist[n].SetBinError(bin, sqrt(binContent));
        }
    }

    // Rescale to the BX
    int config = getConfig();
    switch(config){
        case 9:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-6);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i].GetYaxis()->SetTitle("hits");
        }
        break;

        case 8:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-5);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i].Scale(1E1);
            hitHist[i].GetYaxis()->SetTitle("hits");
        }
        break;

        case 7:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-4);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i].Scale(1E2);
            hitHist[i].GetYaxis()->SetTitle("hits");
        }
        break;

        case 6:
        for(int i=0; i<2; i++){
            edepStrip[2*i].Scale(1E-3);
            edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
            hitHist[i].Scale(1E2);
            hitHist[i].GetYaxis()->SetTitle("hits");
        }
        break;

        default:
        cout << "Config not found. getConfig() = " << config << endl;
    }

    // Fitting
    gStyle->SetOptFit(1);
    TString verbQuiet = "";
    if(verbosityLevel < 0) verbQuiet = "0Q";
    for(Int_t i=0; i<2; i++){
        TF1 fit = TF1("fit", "gaus", fitRange[0], fitRange[1]);
        fit.SetParameter(0, hitHist[i].GetMaximum());
        fit.SetParameter(1, hitHist[i].GetMean());
        fit.SetParameter(2, hitHist[i].GetRMS());
        hitHist[i].Fit("fit", "R"+verbQuiet);
        rchisq[i] = fit.GetChisquare();
    }

    if(verbosityLevel > 0){
        // Output plots
        TCanvas* fitHitsCanvas = new TCanvas("fitHitsCanvas", "Plots", 1000, 600);
        fitHitsCanvas->DivideSquare(4);
        for(Int_t i=0; i<2; i++){
            // fitHitsCanvas->cd(1+i);
            // edepStrip[i*2].Draw("hist");

            fitHitsCanvas->cd(3+i);
            hitHist[i].DrawClone("e");
        }

        // Save histograms to file  
        TString spec = "";
        switch(config){
            case 9:
            spec += "9_";

            case 7:
            if(ifilename.Contains("_sbr")){
                spec = "sbr_";
            }else{
                spec += "7_";
            }
        }
        if(ifilename.Contains("_r12_")) spec += "_r12_";
        // auto foutput = new TFile(spec + "output.root", "UPDATE");
        // edepStrip[0].Write();
        // edepStrip[2].Write();
        // hitHist[0].Write();
        // hitHist[1].Write();
        // foutput->Close();
        // delete foutput;
    }
    // Return chi square
    return rchisq;
}
// Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Default threshold of 1keV and range 90,110.
double* fitHitsChi2(double threshold = 1.){
    int defRange[2] = {90, 110};
    return fitHitsChi2(threshold, defRange);
}
// 2D-graph of the chisquare as a function of threshold and range
std::vector<TGraph2D> fitVSThrRng(){
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1;
    TGraph2D* chiSGraph[2];
    for(int i=0; i<2; i++){
        chiSGraph[i] = new TGraph2D();
        TString tmp = (i ? "Upstream": "Downstream");
        chiSGraph[i]->SetTitle(tmp + " det. #chi^{2}(strip_{window}, e_{thr}) threshold; fit window [#Delta strip]; Threshold [keV]; #chi^{2}");
    }
    for(double tS=0; tS<=10; tS++){
        double threshold = 1.0 + (tS/10)*10.0;
        for(int sw=10; sw<50; sw++){
            int window[2] = {100 - sw, 100 + sw};
            double* chiArr = fitHitsChi2(threshold,window);
            chiSGraph[0]->AddPoint(sw, threshold,chiArr[0]);
            chiSGraph[1]->AddPoint(sw, threshold,chiArr[1]);
        }
    }
    verbosityLevel = verbosityLevel_bak;
    
    std::vector<TGraph2D> retVar;
    for(int i=0; i<2; i++){
        retVar.push_back(*chiSGraph[i]);
        delete chiSGraph[i];
    }
    return retVar;
}
// Plot fitVSThrRng()
void plot_fitVSThrRng(){
    std::vector<TGraph2D> chiSGraph = fitVSThrRng();
    gStyle->SetPalette(1);
    new TCanvas();
    chiSGraph[0].DrawClone("surf1");
    new TCanvas();
    chiSGraph[1].DrawClone("surf1");
}
/*
* ******************************************************************************************************************************************
* Create report with all the major analysis - 1.5/3 minutes
* ******************************************************************************************************************************************
*/
// Save plots in the folder plots
int createReport(TString path = "plots/", TString format = ".pdf", TString ofilename = "out.root"){
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1;

    TFile* ofile = new TFile(ofilename, "recreate");
    TCanvas* painter = new TCanvas("painter", "Title");
    painter->SetGridx(); painter->SetGridy();
    TString fname;

    std::vector<TH1D> temp;
    std::vector<TH2D> temp2;
    // PRIMARY BEAM
    cout << "beamMonitor OK\n";
    temp = beamMonitor();
    for(int i=0; i < 3; i++){
        temp[i].Write();
        temp[i].Draw();
        fname = temp[i].GetName();
        painter->Print(path+fname+format);
    }
    temp.clear();
    
    cout << "primaryProfile OK\n";
    temp2 = primaryProfile();
    temp2[0].Write();
    temp2[0].Draw("LEGO1");
    fname = temp2[0].GetName();
    painter->Print(path+fname+format);
    temp2.clear();
    
    cout << "primaryEnergyProfile OK\n";
    temp2 = primaryEnergyProfile();
    temp2[0].Write();
    temp2[0].Draw("LEGO2Z");
    fname = temp2[0].GetName();
    painter->Print(path+fname+format);
    temp2.clear();
    
    cout << "primarySpectrum OK\n";
    temp = primarySpectrum(11, 16.0);
    for(int i=0; i<2; i++){
        temp[i].Write();
    }
    TLegend* legend = new TLegend(1-0.3,0.7,1-0.1,0.9);
    legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
    legend->AddEntry(&temp[1],TString::Format("#gamma,e^{-},e^{+} - %i entries",(int)temp[1].GetEntries()),"l");
    legend->AddEntry(&temp[0],TString::Format("e^{-},e^{+}    - %i entries",(int)temp[0].GetEntries()),"l");
    temp[0].SetStats(0); temp[1].SetStats(0);
    temp[1].Draw(); temp[0].Draw("same");
    legend->Draw();
    fname = temp[0].GetName();
    painter->Print(path+fname+format);
    painter->SetLogx(); painter->SetLogy();
    painter->Print(path+fname+"LogLog"+format);
    painter->SetLogx(0); painter->SetLogy(0);
    delete legend;
    temp.clear();
    
    cout << "primarySpectrumLogLog OK\n";
    temp = primarySpectrumLogLog(11, 1E-9, 16.0, 100);
    for(int i=0; i<2; i++){
        temp[i].Write();
    }
    TLegend* legend4LogLog = new TLegend(0.1,0.7,0.4,0.9);
    legend4LogLog->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
    legend4LogLog->AddEntry(&temp[1],TString::Format("#gamma,e^{-},e^{+} - %i entries",(int)temp[1].GetEntries()),"l");
    legend4LogLog->AddEntry(&temp[0],TString::Format("e^{-},e^{+}    - %i entries",(int)temp[0].GetEntries()),"l");
    temp[0].SetStats(0); temp[1].SetStats(0);
    painter->SetLogx(); painter->SetLogy();
    temp[1].Draw(); temp[0].Draw("same");
    legend4LogLog->Draw();
    fname = temp[0].GetName();
    painter->Print(path+fname+format);
    painter->SetLogx(0); painter->SetLogy(0);
    delete legend4LogLog;
    temp.clear();


    // Replaced by edepSpectrumChg
    // cout << "edepSpectrum OK\n";
    // temp = edepSpectrum(1000, 0, 1000);
    // for(int i=0; i<2; i++){
    //     temp[i].Write();
    //     temp[i].Draw();
    //     fname = temp[i].GetName();
    //     painter->Print(path+fname+format);
    // }
    // temp.clear();
    
    
    cout << "edepSpectrumChg OK\n";
    temp = edepSpectrumChg(1000, 0, 1000);
    for(int i=0; i < 2; i++){
        TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
        legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
        legend->AddEntry(&temp[i],TString::Format("#gamma,e^{-},e^{+} - %i entries",(int)temp[i].GetEntries()),"l");
        legend->AddEntry(&temp[i+2],TString::Format("e^{-},e^{+}    - %i entries",(int)temp[i+2].GetEntries()),"l");
        temp[i].SetStats(0); temp[i+2].SetStats(0);
        temp[i].Draw();
        temp[i+2].Draw("same");
        legend->Draw();
        fname = temp[i+2].GetName();
        painter->Print(path+fname+format);
        //LogScale version
        legend->SetFillStyle(0);
        painter->SetLogx(1); painter->SetLogy(1);
        painter->Print(path+fname+"LogLog"+format);
        painter->SetLogx(0); painter->SetLogy(0);
        delete legend;
    }
    temp.clear();

    cout << "edepStripChg OK\n";
    temp = edepStripChg();
    for(int i=0; i<2; i++){
        painter->SetName("edepStripChg"+whichDetShort(i));
        painter->SetTitle("Energy deposited per strip in the " + whichDetVerbose(i) + " detector");
        TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
        legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
        legend->AddEntry(&temp[i], TString::Format("#gamma,e^{-},e^{+} - %i entries",(int)temp[i].GetEntries()),"l");
        legend->AddEntry(&temp[i+2], TString::Format("e^{-},e^{+}    - %i entries",(int)temp[i+2].GetEntries()),"l");
        temp[i].SetStats(0); temp[i+2].SetStats(0);
        temp[i].Draw("hist");
        temp[i+2].Draw("hist same");
        legend->Draw();
        fname = painter->GetName();
        painter->Print(path+fname+format);
        delete legend;
    }
    temp.clear();


    cout << "eSpectrumStrip OK\n";
    temp = eSpectrumStrip(/*int stripNb = */100, /*float range = */60.);
    for(int i=0; i<2; i++){
        temp[i].Write();
        temp[i].Draw();
        fname = temp[i].GetName();
        painter->Print(path+fname+format);
    }
    temp.clear();

    cout << "eLongitudinal_array OK\n";
    temp = eLongitudinal_array(/*nbBins = */100);
    for(int i=0; i < 2; i++){
        temp[i].Write();
        temp[i].Draw("hist");
        fname = temp[i].GetName();
        painter->Print(path+fname+format);
    }
    temp.clear();

    cout << "edepProject OK\n";
    temp = edepProject(/*meshSize = */0.100);
    for(int i=0; i < 4; i++){
        painter->cd(i);
        temp[i].Write();
        temp[i].Draw("hist");
        fname = temp[i].GetName();
        painter->Print(path+fname+format);
    }
    temp.clear();

    cout << "doseXY OK\n";
    temp2 = doseXY();
    for(int i=0; i < 4; i++){
        temp2[i].Write();
        temp2[i].Draw("colz");
        fname = temp2[i].GetName();
        painter->Print(path+fname+format);
    }
    temp2.clear();

    cout << "fitHits OK\n";
    std::vector<TH1I> fitHitsPlots = fitHits();
    for(int i=0; i < 2; i++){
        fitHitsPlots[i].Write();
        fitHitsPlots[i].Draw();
        fname = fitHitsPlots[i].GetName();
        painter->Print(path+fname+format);
    }
    fitHitsPlots.clear();
    
    cout << "edepProjectPrimEnPartitioned OK\n";
    // Lambda function for formatting
    auto formatEnergy = [](double x){
        TString unit;
        if(x>=0.1){
            unit = TString::Format("%.1f GeV", x);
        }else if(0.001<=x && x<0.1){
            unit = TString::Format("%.1f MeV", x*1E3);
        }else if(0.000001<=x && x<0.001){
            unit = TString::Format("%.1f keV", x*1E6);
        }else if(x<0.000001){
            unit = TString::Format("%.1f eV", x*1E9);
        }else{
            unit = TString::Format("%f", x);
        }
        return unit;
    };
    std::vector<TH1D> partitions = edepProjectPrimEnPartitioned();
    int slices = partitions.size() / 2 - 1; //size is always even, so there is no possibility for int rounding errors
    for(int det=0; det<2; det++){
        double maxCounts = 0;
        TLegend* edepProjectPrimEnPartitionedLegend = new TLegend(1.0-0.4, 0.7, 1.0-0.1, 0.9);
        edepProjectPrimEnPartitionedLegend->AddEntry(&partitions[2*slices+det], partitions[2*slices+det].GetTitle(), "l");//"0 < E < " + formatEnergy(partitions[2*slices+det].GetMaximum())
        for(int i= slices-1; i>=0; i--){
            edepProjectPrimEnPartitionedLegend->AddEntry(&partitions[2*i+det], partitions[2*i+det].GetTitle(), "l");
        }
        // Make sure that all histograms are viewed correctly
        for(int i= slices-1; i>=0; i--){
            if(partitions[2*i+det].GetMaximum() > maxCounts) maxCounts = partitions[2*i+det].GetMaximum();
        }
        partitions[2*slices+det].SetStats(0);
        partitions[2*slices+det].SetMinimum(1);
        partitions[2*slices+det].Draw("hist");
        for(int i= slices-1; i>=0; i--){
            partitions[2*i+det].GetYaxis()->SetRangeUser(0, 1.05*maxCounts); // Make sure that all histograms are viewed correctly
            partitions[2*i+det].SetStats(0);
            partitions[2*i+det].SetMinimum(1);
            partitions[2*i+det].Draw("hist same");
        }

        edepProjectPrimEnPartitionedLegend->Draw();
        fname = "edepProjectPrimEnPartitioned"+whichDetShort(det);
        painter->Print(path+fname+format);
    }
    partitions.clear();

    // cout << "fitVSThrRng OK\n";
    // std::vector<TGraph2D> fitVSThrRngPlots = fitVSThrRng();
    // for(int i=0; i < 2; i++){
    //     fitVSThrRngPlots[i].Write();
    //     gStyle->SetPalette(1);
    //     fitVSThrRngPlots[i].Draw("surf1");
    //     fname = fitVSThrRngPlots[i].GetName();
    //     painter->Print(path+fname+format);
    // }
    // fitVSThrRngPlots.clear();

    ofile->Close();
    verbosityLevel = verbosityLevel_bak;
    return 0;                                              
}


/*
* ******************************************************************************************************************************************
* NEW METHODS!!!
* ******************************************************************************************************************************************
*/

// Superimpose the beam profile of the primaries with the reconstructed one
void overlayReconstructedBeamProfile(){
    TH2D primaryEPro = primaryEnergyProfile()[0];
    TH1D* primProjX = primaryEPro.ProjectionX();

    double scale;
    scale = primProjX->GetXaxis()->GetBinWidth(1)/(primProjX->Integral());
    primProjX->Scale(scale);
    primProjX->Draw("hist");

    std::vector<TH1D> detEPro = edepProject();
    TH1D det0XProj = detEPro[0];
    scale = det0XProj.GetXaxis()->GetBinWidth(1)/(det0XProj.Integral());
    det0XProj.Scale(scale);
    det0XProj.SetLineColor(kRed);
    det0XProj.DrawClone("hist same");
}
/*
* ******************************************************************************************************************************************
* Count events around an axis, assuming radial symmetry - START
* ******************************************************************************************************************************************
*/
// This function explains the choose for the radius used in the average dose calculation.
void toolsCumulativeProbabilityXYGaussian(double events = 1E5, float sig = 0.375, float mult = 1.){
    TH2D* h = new TH2D("h", "h", 100, -2, 2, 100, -2, 2);
    TRandom3* rand = new TRandom3(3);

    TGraph* sigmaPlot = new TGraph();
    sigmaPlot->SetTitle(TString::Format("Cumulative probability - sigma %f - events %f", sig, events));
    sigmaPlot->GetXaxis()->SetTitle("\\frac{r}{\\sigma}");
    sigmaPlot->GetYaxis()->SetTitle("Probability");

    double ePoint = mult + 5.;
    TGraph* sigmaFit = new TGraph();

    for(double sPoint = mult; sPoint < ePoint; sPoint+=0.2){
        double counter = 0;
        double radius2 = 2. * pow( sPoint * sig, 2.);
        // double radius = sig * sqrt( 2.*log(20.));
        // radius2 = radius*radius;
        for(double i=0; i< events; i++){
            double x = rand->Gaus(0, sig);
            double y = rand->Gaus(0, sig);
            // y=0;
            // h->Fill(x,y);
            if(x*x + y*y < radius2) counter++;
        }
        double sigValue = counter/events;
        sigmaPlot->AddPoint(sqrt(radius2)/sig, sigValue);
        sigmaFit->AddPoint(sqrt(radius2)/sig, 1 - exp(-radius2/(2.*sig*sig)) );
        // cout << "counter/total: " << value << endl;
    }
    sigmaPlot->SetLineWidth(0);
    sigmaPlot->Draw("AP*");
    sigmaFit->Draw("same");
}
void macro(){}

/*
* ******************************************************************************************************************************************
* Utilities
* ******************************************************************************************************************************************
*/
// void exportBruschi(){
//     TCanvas* resultCanvas = new TCanvas("resultCanvas", "Plots", 1000, 600);
//     resultCanvas->Divide(2,2);
    
//     loadAnotherFile("kyle/kyleBeam_9_50_0_kapton.root");

//     doseXY(0, 0.100);
//     edepXYHist->Scale(1E-6);

//     resultCanvas->cd(1);
//     auto projXDet0 = edepXYHist->ProjectionX();
//     projXDet0->SetTitle("det0X");
//     projXDet0->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
//     projXDet0->DrawClone("hist");
//     exportDAT(projXDet0, "det0_edepX.txt");
//     exportCSV(projXDet0, "det0_edepX.csv");

//     resultCanvas->cd(2);
//     auto projYDet0 = edepXYHist->ProjectionY();
//     projYDet0->SetTitle("det0Y");
//     projYDet0->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
//     projYDet0->DrawClone("hist");
//     exportDAT(projYDet0, "det0_edepY.txt");
//     exportCSV(projYDet0, "det0_edepY.csv");
    

//     doseXY(1, 0.100);
//     edepXYHist->Scale(1E-6);

//     resultCanvas->cd(3);
//     auto projXDet1 = edepXYHist->ProjectionX();
//     projXDet1->SetTitle("det1X");
//     projXDet1->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
//     projXDet1->DrawClone("hist");
//     exportDAT(projXDet1, "det1_edepX.txt");
//     exportCSV(projXDet1, "det1_edepX.csv");

//     resultCanvas->cd(4);
//     auto projYDet1 = edepXYHist->ProjectionY();
//     projYDet1->SetTitle("det1Y");
//     projYDet1->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
//     projYDet1->DrawClone("hist");
//     exportDAT(projYDet1, "det1_edepY.txt");
//     exportCSV(projYDet1, "det1_edepY.csv");
// }
/*
* ******************************************************************************************************************************************
* Old methods - obsolete/unuseful (kept here because they have useful things accesible by the search tool)
* ******************************************************************************************************************************************
*/
// Energy spectrum of depositions for the simulation in air and vacuum at comparison. It uses kyleBeam_9_100 
// void eSpectrumAirVacuum(TString sel = "5"){
//     loadAnotherFile("kyle/kyleBeam_6_"+sel+"0_0.root");
//     TH1D air[2], vacuum[2];
//     for(int i=0; i<2; i++){
//         air[i] = edepSpectrum(i)[0];
//         air[i].SetStats(0);
//     }

//     loadAnotherFile("kyle/kyleBeam_6_"+sel+"0_0_vacuum.root");
//     for(int i=0; i<2; i++){
//         vacuum[i] = edepSpectrum(i)[0];
//         vacuum[i].SetLineColor(kRed);
//         vacuum[i].SetStats(0);
//     }

//     TLegend* legend[2];
//     for(int i=0; i<2; i++){
//         TCanvas* compCanvas = new TCanvas("eSpectrumAirVacuum"+whichDetShort(i), "Comparison");
//         vacuum[i].DrawClone();
//         air[i].DrawClone("same");

//         legend[i] = new TLegend(1.0-0.4, 0.7, 1.0-0.1, 0.9);
//         legend[i]->AddEntry(&air[i], "air", "l");
//         legend[i]->AddEntry(&vacuum[i], "vacuum", "l");
//         legend[i]->DrawClone();
//     }
// }
// Impact of charged particles to the spectrum of deposited energy. Return a 2-array with the spectrum of energy deposited in upstream/downstream detector in the range [a,b] with #bins. Sel is a selector for the values \xi={"5.0", "7.0", "10.0"} 
// int eSpectrumChgCanvas(int binsIn, double a, double b, TString path = "plots/", TString format = ".pdf"){
//     static const int bins = binsIn;
//     TH1D eSpectrChgHist[2];
//     for (int i=0; i<2; i++){
//         eSpectrChgHist[i] = TH1D(TString::Format("eSpectrChgDet%i",i), "Energy spectrum " + whichDetVerbose(i) + " detector", bins, a, b);
//     }

//     // Get the list of tracks of primaries with energy in the given range
//     std::vector<int> pChgList;
//     int eventSel, pdg;
//     primaryTree->SetBranchAddress("event", &eventSel);
//     primaryTree->SetBranchAddress("pdg", &pdg);
//     Long64_t primaryEntries = primaryTree->GetEntries();
//     if(verbosityLevel>0) cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
//     resetProgress();
//     pChgList.reserve(primaryEntries);
//     for(Long64_t i=0; i< primaryEntries; i++ ){
//         float status = (float)(i+1) / primaryEntries;
//         printProgress(status);
        
//         primaryTree->GetEntry(i);
//         if(pdg != 22) pChgList.push_back(i);
//     }

//     Long64_t chgNb = pChgList.size();
//     if(verbosityLevel>0) cout << "Charged primaries: " << chgNb << endl;
//     // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
//     std::sort(pChgList.begin(), pChgList.end());
//     // Select energy deposition events with charged primaries
//     int event, det; double edep;
//     eventTree->SetBranchAddress("event", &event);
//     eventTree->SetBranchAddress("det", &det);
//     eventTree->SetBranchAddress("edep", &edep);
//     Long64_t eventEntries = eventTree->GetEntries();
//     if(verbosityLevel>0) cout << "eventTree has " << eventEntries << " entries. Looping over entries...\n";
//     // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
//     int p0 = 0;
//     resetProgress();
//     for(Long64_t i=0; i< eventEntries; i++ ){
//         eventTree->GetEntry(i);
//         float status = (float)(i+1) / eventEntries;
//         printProgress(status);
//         for(int i=p0; i<chgNb; i++){
//             double val=pChgList[i];
//             if(event < val){
//                 break;
//             }else if(event == val){
//                 eSpectrChgHist[det].Fill(edep);
//                 p0=i;
//                 break;
//             }
//         }
//     }
//     // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;

//     // cout << "Entry point. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
//     int verbosityLevel_bak = verbosityLevel;
//     verbosityLevel = -1;
//     std::vector<TH1D> eSpectrTotHist = edepSpectrum(bins, a, b);
//     // cout << "Breakpoint A. Bins, a, b " << bins << '\t' << a  << '\t' << b << endl;
//     // cout << "Breakpoint A1" << endl;
//     verbosityLevel = verbosityLevel_bak;
//     // cout << "Breakpoint B" << endl;
//     for(int i=0; i<2; i++){
//         cout << "Breakpoint C. Loop " << i << endl;
//         eSpectrChgHist[i].GetXaxis()->SetTitle("Energy [keV]");
//         eSpectrChgHist[i].GetYaxis()->SetTitle("counts");

//         eSpectrTotHist[i].GetXaxis()->SetTitle("Energy [keV]");
//         eSpectrTotHist[i].GetYaxis()->SetTitle("counts");

//         eSpectrTotHist[i].SetStats(0);
//         eSpectrChgHist[i].SetStats(0);
//         eSpectrChgHist[i].SetLineColor(kRed);
        
//         TCanvas* eSpectrumChgCanvas = new TCanvas("eSpectrumChg"+whichDetShort(i)+"Canvas", "Impact of charged particles to the spectrum of deposited energy");
//         // cout << "Breakpoint D. Loop " << i << endl;
//         TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
//         legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
//         legend->AddEntry(&eSpectrTotHist[i],TString::Format("#gamma,e^{-},e^{+} - %i entries",(int)eSpectrTotHist[i].Integral()),"l");
//         legend->AddEntry(&eSpectrChgHist[i],TString::Format("e^{-},e^{+}    - %i entries",(int)eSpectrChgHist[i].Integral()),"l");
//         eSpectrTotHist[i].Draw("hist");
//         eSpectrChgHist[i].Draw("same hist");
//         // cout << "Breakpoint E. Loop " << i << endl;
//         legend->Draw();
//         cout << "Breakpoint F. Loop " << i << endl;
//         TString fname = "eSpectrumChg"+whichDetShort(i);
//         eSpectrumChgCanvas->Print(path+fname+format);
//         delete legend, eSpectrumChgCanvas;
//         cout << "Breakpoint C. End loop " << i << endl;
//     }
//     pChgList.clear();
//     cout << "Breakpoint G" << endl;
//     return 0;
// }
// int superimpose(double thr = 1.0){
//     loadAnotherFile("kyle/kyleBeam_7_50_0_kapton.root");
//     detectorPerformanceBX(thr);
//     TH1D* primaryHistX = new TH1D("primaryHistX", "Edep/BX upstream detector;strip number; (normalized)", 200, 1, 200);
//     TH1D* primaryHistY = new TH1D("primaryHistY", "Edep/BX downstream detector;strip number; (normalized)", 200, 1, 200);
//     // Initialize variables for the calculation
//     double x0, y0, ekin;
//     cout << "Setting branch addresses...";
//     primaryTree->SetBranchAddress("x0", &x0);                       // mm
//     cout << "x0,";
//     primaryTree->SetBranchAddress("y0", &y0);   
//     cout << "y0,";
//     primaryTree->SetBranchAddress("ekin", &ekin);   
//     cout << "ekin.\n";

//     Long64_t primaryEntries = primaryTree->GetEntries();
//     cout << "primaryTree has " << primaryEntries << " entries. Looping over entries...\n";
//     resetProgress();
//     for(Long64_t i=0; i< primaryEntries; i++ ){
//         float status = (float)(i+1) / primaryEntries;
//         printProgress(status);
//         primaryTree->GetEntry(i);

//         int stripXNb = (x0 + 10.0) / 20. * 200. + 1.;
//         primaryHistX->Fill(stripXNb, ekin);                                    // Fill the total energy deposition histogram
//         int stripYNb = (y0 + 10.0) / 20. * 200. + 1.;
//         primaryHistY->Fill(stripYNb, ekin);                                    // Fill the total energy deposition histogram
//     }
    

//     new TCanvas();

//     double scale;
//     scale = primaryHistX->GetXaxis()->GetBinWidth(1)/(primaryHistX->Integral());
//     primaryHistX->Scale(scale);
//     primaryHistX->SetLineColor(kRed);
//     primaryHistX->Draw("hist");

//     scale = energyHist[0]->GetXaxis()->GetBinWidth(1)/(energyHist[0]->Integral());
//     energyHist[0]->Scale(scale);
//     energyHist[0]->Draw("hist same");

//     // new TCanvas();
//     // primaryHistY->Draw("hist");

//     return 0;
//     // new TCanvas();
//     // energyHist[0]->Draw("hist");
//     // primaryHist[0]->SetLineColor(kRed);
//     // primaryHist[0]->Draw("same hist");
    

//     // new TCanvas();
//     // energyHistCopy[1]->SetLineColor(kRed);
//     // energyHistCopy[1]->Draw();
//     // energyHist[1]->Draw("same");
// }
// TH1I* hitHist[2];
// TH1D* energyHist[2];
// int estimateDetectorPerformance(double threshold = 1.){
//     // Prepare the histogram
//     // 1 bin per strip from strip 1 to 200
//     energyHist[0] = new TH1D("energyX", "Edep upstream detector;strip number; energy [keV];", 200, 1, 200);
//     energyHist[0]->GetXaxis()->SetRangeUser(60,140);
//     hitHist[0] = new TH1I("hitX", "Hits in the upstream detector;Strip number [1-200];counts", 200, 1, 200);
//     hitHist[0]->GetXaxis()->SetRangeUser(60,140);

//     energyHist[1] = new TH1D("energyY", "Edep downstream detector;strip number; energy [keV];", 200, 1, 200);
//     energyHist[1]->GetXaxis()->SetRangeUser(60,140);
//     hitHist[1] = new TH1I("hitY", "Hits in the downstream detector;Strip number [1-200];counts", 200, 1, 200);
//     hitHist[1]->GetXaxis()->SetRangeUser(60,140);

//     // Initialize variables for the calculation
//     int det, strip;
//     double edep;
//     cout << "Setting branch addresses...";
//     stripTree->SetBranchAddress("det", &det);                  // 0,1
//     cout << "det,";
//     stripTree->SetBranchAddress("strip", &strip);              // 0-200
//     cout << "strip,";
//     stripTree->SetBranchAddress("edep", &edep);
//     cout << "edep.\n";

//     double totEnergy[2] = {0};
//     Long64_t stripEntries = stripTree->GetEntries();
//     cout << "stripTree has " << stripEntries << " entries. Looping over entries...\n";
//     resetProgress();
//     for(Long64_t i=0; i< stripEntries; i++ ){
//         float status = (float)(i+1) / stripEntries;
//         printProgress(status);
//         stripTree->GetEntry(i);

//         totEnergy[det] += edep;                                                // Total energy, useful for normalization (optional)
//         energyHist[det]->Fill(strip, edep);                                    // Fill the total energy deposition histogram
        
//         if(edep < threshold) continue;                                         // Hit condition. If edep is less than threshold keV then quit the loop istance
        
        
//         hitHist[det]->Fill(strip, (int)(edep/threshold));                      // Add the number of hits in the corresponding x/y det hist.
//                                                                                // For example: edep=5.4 keV -> 5 hits
//     }

//     // I've noticed that error are not properly calculated, despite the option is
//     //      root [13] result->GetBinErrorOption()
//     //      (TH1::EBinErrorOpt) (TH1::kNormal) : (int) 0)
//     // correct. So I wrote the following code to fix the issue
//     cout << "Fixing errors...\n";
//     resetProgress();
//     for(int n=0; n < 2 ; n++ ){
//         int binNb = hitHist[n]->GetNbinsX();
//         for(int bin=1; bin<=binNb; bin++){
//             float status = (float)(bin*(n+1)-1)/(float)(2 * binNb);
//             printProgress(status);

//             double binContent = hitHist[n]->GetBinContent(bin);
//             hitHist[n]->SetBinError(bin, sqrt(binContent));
//         }
//     }


//     TCanvas* c1 = new TCanvas("c1", "Plots", 1000, 600);
//     c1->Divide(2,2);

//     c1->cd(1);
//     energyHist[0]->Draw("hist");

//     c1->cd(2);
//     energyHist[1]->Draw("hist");
    
    
//     gStyle->SetOptFit(1);
//     TF1* fit = new TF1("fit", "gaus", 90, 110);
//     TF1* fit2 = new TF1("fit2", "gaus", 90, 110);

//     c1->cd(3);
//     hitHist[0]->Draw("e");
//     fit->SetParameter(0, hitHist[0]->GetMaximum());
//     fit->SetParameter(1, hitHist[0]->GetMean());
//     fit->SetParameter(2, hitHist[0]->GetRMS());
//     hitHist[0]->Fit("fit", "R");

//     c1->cd(4);
//     hitHist[1]->Draw("e");
//     fit2->SetParameter(0, hitHist[1]->GetMaximum());
//     fit2->SetParameter(1, hitHist[1]->GetMean());
//     fit2->SetParameter(2, hitHist[1]->GetRMS());
//     hitHist[1]->Fit("fit2", "R");
    

//     TString spec = "";
//     if(ifilename.Contains("_sbr")){
//         spec = "sbr_";
//     }
//     if(ifilename.Contains("_7_")){
//         spec += "7_";
//     }else if(ifilename.Contains("_9_")){
//         spec += "9_";
//     }
//     if(ifilename.Contains("_r12_")) spec += "_r12_";
//     // Save histograms to file
//     TFile* output = new TFile(spec + "output.root", "UPDATE");
//     energyHist[0]->Write();
//     energyHist[1]->Write();
//     hitHist[0]->Write();
//     hitHist[1]->Write();
//     output->Close();
//     delete output;

//     // Save image of the canvas
//     c1->Print(spec + "detectorPerformances.png");

//     if(verbosityLevel>0){
//         cout << "-----------------------------------------------------------------" << endl;
//         cout << "File: " << ifilename << endl;
//         cout << "-----------------------------------------------------------------" << endl;
//     }
//     // Dose calculations
//     for(int i=0; i<2; i++){
//         double conversionFactor = 0.1 * 1.60 / 3.970;
//         double peakEdepStrip = energyHist[i]->GetMaximum();
//         double peakStrip = energyHist[i]->GetMaximumBin();
//         double vol = 100 * 100 * 100; //um3
//         double peakDoseStrip = peakEdepStrip / vol * conversionFactor;
//         if(verbosityLevel>1){
//             cout    << "Detector: "                         << whichDetVerbose(i) << endl;
//             if(ifilename.Contains("_7_")){
//                 cout    << "Total energy deposit is:\t"             << totEnergy[i]/1E6 << " GeV" << endl;
//                 cout    << "Total charge is:\t\t\t"                 << totEnergy[i]*1E3 / 27 * 1.6E-19 * 1E12 << " pC (#e/h x e, with #pairs given by: " << totEnergy[i]*1E3 / 27 << " )" << endl;
//                 cout    << "Highest energy dep is:\t\t"             << peakEdepStrip << " keV in strip #" << peakStrip << endl;
//                 cout    << "Peak dose is:\t\t\t"                    << peakDoseStrip << " Gy in strip " << peakStrip << " (energy: " << peakEdepStrip << " keV, vol: " << vol << " um3)" <<  endl;
//                 cout    << "----------------BX--------------------" << endl;
//                 cout    << "Total energy deposit/BX is:\t"          << totEnergy[i]/1E6 * 1E2 << " GeV" << endl;
//                 cout    << "Total charge/BX is:\t\t"                << (totEnergy[i]*1E2)*1E3 / 27 * 1.6E-19 * 1E12 << " pC (#e/h x e, with #pairs given by: " << (totEnergy[i]*1E2)*1E3 / 27 << " )" << endl;
//                 cout    << "Highest energy dep/BX is:\t"            << peakEdepStrip*1E2 << " keV in strip #" << peakStrip << endl;
//                 cout    << "Peak dose/BX is:\t\t"                   << peakDoseStrip*1E2 << " Gy in strip " << peakStrip << " (energy: " << peakEdepStrip*1E2 << " keV, vol: " << vol << " um3)" <<  endl;
//             }else if(ifilename.Contains("_9_")){
//                 cout    << "Total energy deposit/BX is:\t\t"        << totEnergy[i]/1E6 << " GeV" << endl;
//                 cout    << "Total charge/BX is:\t\t\t"              << totEnergy[i]*1E3 / 27 * 1.6E-19 * 1E12 << " pC (#e/h x e, with #pairs given by: " << totEnergy[i]*1E3 / 27 << " )" << endl;
//                 cout    << "Highest energy dep/BX is:\t\t"          << peakEdepStrip << " keV in strip #" << peakStrip << endl;
//                 cout    << "Peak dose/BX is:\t\t\t"                 << peakDoseStrip << " Gy in strip " << peakStrip << " (energy: " << peakEdepStrip << " keV, vol: " << vol << " um3)" <<  endl;
//             }
//             cout    << "*****************************************************************************\n\n";
//             }
//     }

//     if(ifilename.Contains("_7_")){
//         return 7;
//     }else if(ifilename.Contains("_9_")){
//         return 9;
//     }else{
//         return -1;
//     }
// }
// // Print energy deposited per strip rinormalized to BX. Vertical energy unit is GeV
// void detectorPerformanceBX(double threshold = 1.){
//     TString spec = whichConfigIs();
//     int result = estimateDetectorPerformance(threshold);

//     if(result == -1){
//         cout << "WARNING: an error as occurred!" << endl;
//         return;
//     }else if(result == 7){
//         energyHist[0]->Scale(1E2);
//         energyHist[1]->Scale(1E2);
//     }

//     energyHist[0]->Scale(1E-6);     //keV -> GeV
//     energyHist[1]->Scale(1E-6);

//     energyHist[0]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");    //Set right axis name
//     energyHist[1]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
//     energyHist[0]->SetTitle("Edep in strips upstream detector");
//     energyHist[1]->SetTitle("Edep in strips downstream detector");

//     TCanvas* edepStripBXUpCanvas = new TCanvas("edepStripBXUp", "edepStripBXUp");
//     energyHist[0]->Draw("hist");
//     edepStripBXUpCanvas->Print("det0_" + spec + "edepBXStrip.png");

//     TCanvas* edepStripBXDoCanvas = new TCanvas("edepStripBXDo", "edepStripBXDo");
//     energyHist[1]->Draw("hist");
//     edepStripBXDoCanvas->Print("det1_" + spec + "edepBXStrip.png");

//     TCanvas* c1 = new TCanvas("c1", "Plots", 1000, 600);
//     c1->Divide(2,2);
//     c1->cd(1);
//     energyHist[0]->Draw("hist");
//     c1->cd(2);
//     energyHist[1]->Draw("hist");
//     c1->cd(3);
//     hitHist[0]->Draw("e");
//     c1->cd(4);
//     hitHist[1]->Draw("e");
//     c1->Print(spec + "detectorperformanceBX.png" );

//     // Save histograms to file
//     TFile* output = new TFile(spec + "BXoutput.root", "UPDATE");
//     energyHist[0]->Write();
//     energyHist[1]->Write();
//     hitHist[0]->Write();
//     hitHist[1]->Write();
//     output->Close();
//     delete output;
// }
// Compare total energy depositions with the standalone MC 1D
// void compareTotEdep(int detNb=0, int bins = 250, float emax = 500){
//     TString plotTitle = "Energy spectrum in the " + whichDetVerbose(detNb) + " detector";
//     TH1F* totEdep = new TH1F("totEdep", plotTitle + ";Deposited energy [keV];", bins, 0, emax);
//     eventTree->Project("totEdep", "edep", TString::Format("det==%i",detNb) );    
//     totEdep->SetLineColor(kRed);
//     totEdep->SetStats(0);

//     loadAnotherFile("data/run_sbr.root");
//     TH1F* totEdepSbr = new TH1F("totEdepSbr", plotTitle + ";Deposited energy [keV];", bins, 0, emax);
//     eventTree->Project("totEdepSbr", "edep", TString::Format("det==%i",detNb) );
//     totEdepSbr->SetLineColor(kBlue);
//     totEdepSbr->SetStats(0);

//     TCanvas* canvas = new TCanvas("plotEdepDet", plotTitle);
//     totEdep->Draw();
//     totEdepSbr->Draw("same");
//     canvas->Print(whichDetShort(detNb) + "_edepComparison.png");
// }

// For tiny mesh, the array size greatly exceeds the memory capability
// therefore the "coarse graining" process in divided in N jobs
// For each job, a mesh grid is applied and peak dose is found. Then
// the resulting maxima are compared to estimate the highest peak dose.
// double doseForTinyMesh(int detNb=0, double meshSize = 0.100 /*mm*/, double range = 20. /*mm*/){
//     int memSize = 1000;
//     double nbMeshes = range/meshSize;
//     int nbJobs = nbMeshes / memSize + 1; //adds a job for making sure that every data is taken
    
//     double conversionFactor = 0.1 * 1.60 / 3.970;
//     double vol = (meshSize*1E3)*(meshSize*1E3)*(100); //um3
//     double peakDose = 0;
//     double jobPeakDose = 0;
    
//     double chunkSize = range / nbMeshes * memSize;
//     double chunkX = -range/2. + chunkSize/2.;     // mm
//     double chunkY = -range/2. + chunkSize/2.;     // mm
//     cout << "Mesh size: " << meshSize*1E3 << " um | Chunk size: " << chunkSize << " mm | Jobs: " << nbJobs << endl;
    
//     int stats = 1;
//     for(int jobY=0; jobY < nbJobs; jobY++ ){
//         for(int jobX=0; jobX < nbJobs; jobX++ ){
//             // Debug
//             cout    << "Job: " << stats << "/" << nbJobs*nbJobs;

//             int n=memSize;
//             int m=memSize;
//             // if(jobX == nbJobs-1){
//             //     n = chunkSize;
//             // }else if(jobY == nbJobs-1){
//             //     m = memSize;
//             // }
//             vector<vector<double>> samplingBlock;
//             samplingBlock.resize(n,vector<double>(m));
//             for(int i=0;i<n;i++){
//                 for(int j=0;j<m;j++){
//                         samplingBlock[i][j] = 0.;
//                     }
//             }

//             double xmin = chunkX - chunkSize/2.;            double xmax = chunkX + chunkSize/2.;
//             double ymin = chunkY - chunkSize/2.;            double ymax = chunkY + chunkSize/2.;
//             if(xmax > range/2) xmax = range/2; if(ymax > range/2) ymax = range/2;
//             cout << " " << xmin << " " << xmax << " " << ymin << " " << ymax << endl;;

//             int det; double x, y; double edep;
//             debugTree->SetBranchAddress("det", &det);               // 0,1
//             debugTree->SetBranchAddress("x", &x);                   // mm
//             debugTree->SetBranchAddress("y", &y);                   // mm
//             debugTree->SetBranchAddress("edep", &edep);             // keV
            
            
//             cout << "Coarse graining..." << endl;
//             int treeEntries = debugTree->GetEntries();
//             resetProgress();
//             for(int i=0; i< treeEntries; i++ ){
//                 float status = (float)(i+1)/(float)(treeEntries);
//                 printProgress(status);

//                 debugTree->GetEntry(i);
//                 if(det != detNb) continue;                                    // Select the correct detector
//                 // if(x<-1.5 || x>1.5 || y<-1.5 || y>1.5) continue;           // Focus on the subset around the center
//                 if(!(xmin < x && x < xmax && ymin < y && y < ymax)) continue; // Job filtering

//                 int xPosGrained = (x-xmin)/chunkSize * n;                     // find the mesh coordinate in the grid
//                 int yPosGrained = (y-ymin)/chunkSize * m;                     // z
//                 if(xPosGrained == n){                                         // Handle the boundary cases
//                     xPosGrained--;
//                 }else if(yPosGrained == m){
//                     yPosGrained--;
//                 }
//                 // if(xPosGrained < 0 || xPosGrained > n || yPosGrained < 0 || yPosGrained > n) cout << "Something odd." << TString::Format("%i %i", xPosGrained, yPosGrained);
//                 samplingBlock[xPosGrained][yPosGrained] += edep / vol * conversionFactor;               // Fill the corse-grained dose map
//             }

//             // Sampling the peak value
//             cout << "Sampling for the peak value..." << endl;
//             jobPeakDose = 0;
//             resetProgress();
//             for(int j=0; j<m;j++){
//                 for(int i=0; i<n; i++){
//                     float status = (float)(i+j+1)/(float)(n+m);
//                     printProgress(status);

//                     jobPeakDose = max(samplingBlock[i][j], jobPeakDose);
//                 }
//             }

//             cout << TString::Format("peak dose: %f Gy)", jobPeakDose) << endl;
//             peakDose = max(peakDose, jobPeakDose);

//             chunkX += chunkSize;
//             stats++;
//         }
//         chunkX = -range/2 + chunkSize/2.;
//         chunkY += chunkSize;
//     }
//     cout << "Peak dose value: " << peakDose << " Gy"<< endl;
//     return peakDose;
// }
// Dose evaluation.
// double peakDose(int detNb=0, double meshSize = 0.100 /*mm*/){
//     int nbBins = TMath::Nint(20./meshSize);
//     TString plotTitleDose = "Dose XY map in the " + whichDetVerbose(detNb) + " detector";    
//     TH2F* doseXYHist_local = new TH2F("doseXYHist_local", plotTitleDose + " [Gy] ;X [mm]; Y [mm];", nbBins, -10, 10, nbBins, -10, 10);

//     // Debug
//     // cout << "Mesh size [um]: "  << meshSize*1000 << endl;
//     // cout << "Mesh blocks: "     << nbXpts<<'x'<<nbYpts      << endl;
//     // cout << "Hist. bins: "      << doseXYHist_local->GetNbinsX()<<'x'<<doseXYHist_local->GetNbinsY()      << endl;

//     int det, pdg; double x, y; double edep, slen;
//     debugTree->SetBranchAddress("det", &det);               // 0,1
//     debugTree->SetBranchAddress("x", &x);                   // mm
//     debugTree->SetBranchAddress("y", &y);                   // mm
//     debugTree->SetBranchAddress("edep", &edep);             // keV
//     debugTree->SetBranchAddress("pdg", &pdg);               // -11, 11, 22, -

//     double totEnergy = 0;
//     double conversionFactor = 0.1 * 1.60 / 3.970;
//     double vol = (meshSize*1E3)*(meshSize*1E3)*(100);

//     int treeEntries = debugTree->GetEntries();
//     for(int i=0; i< treeEntries; i++ ){
//         debugTree->GetEntry(i);
//         if(det != detNb) continue;                                    // Select the correct detector

//         totEnergy += edep;                                            // Accumulate for the total energy deposition
//         doseXYHist_local->Fill(x, y, edep / vol * conversionFactor);
//     }

//     // Debug (cross-check)
//     Int_t maxBinx, maxBiny, maxBinz;
//     doseXYHist_local->GetBinXYZ(doseXYHist_local->GetMaximumBin(), maxBinx, maxBiny, maxBinz);
//     double maxPosx = doseXYHist_local->GetXaxis()->GetBinCenter(maxBinx);
//     double maxPosy = doseXYHist_local->GetYaxis()->GetBinCenter(maxBiny);
//     double peakDose = doseXYHist_local->GetMaximum();
//     delete doseXYHist_local;

//     return peakDose;
// }

/*
* ******************************************************************************************************************************************
* Sensitivity - START
* ******************************************************************************************************************************************
*/
// void eDepStripVsCut(int detNb=0, int stipNb = 100){
//     TGraph* eDepStripVsCutPlot = new TGraph();
//     eDepStripVsCutPlot->SetTitle(TString::Format("Energy deposited in strip %i - det%i", stipNb, detNb));
//     eDepStripVsCutPlot->GetXaxis()->SetTitle("Cut value [um]");
//     eDepStripVsCutPlot->GetYaxis()->SetTitle("Tot. Energy [keV]");

//     double result =0;
//     loadAnotherFile("data/runXY_7_cut001_air10.root");
//     result = eSpectrumStrip(detNb, stipNb);
//     eDepStripVsCutPlot->AddPoint(0.01, result);

//     loadAnotherFile("data/runXY_7_cut01_air10.root");
//     result = eSpectrumStrip(detNb, stipNb);
//     eDepStripVsCutPlot->AddPoint(0.1, result);

//     loadAnotherFile("data/runXY_7_cut1_air10.root");
//     result = eSpectrumStrip(detNb, stipNb);
//     eDepStripVsCutPlot->AddPoint(1.0, result);

//     loadAnotherFile("data/runXY_7_cut100_air10.root");
//     result = eSpectrumStrip(detNb, stipNb);
//     eDepStripVsCutPlot->AddPoint(100, result);

//     loadAnotherFile("data/run_sbr.root");
//     result = eSpectrumStrip(detNb, 95);
//     eDepStripVsCutPlot->AddPoint(1000, result);

//     eDepStripVsCutPlot->Draw();
// }

// Return the graph for the study of the peak dose sensitivity
// TGraph* peakDoseMeshSensitivity(int detNb=0, int points = 100, float maxMesh = 1.){
//     TGraph* pkSensiPlot = new TGraph();
//     pkSensiPlot->SetTitle("Peak dose sensitivity respect to mesh size");
//     pkSensiPlot->GetXaxis()->SetTitle("Mesh size L [um]");
//     pkSensiPlot->GetYaxis()->SetTitle("Peak dose [mGy]");
    
//     // Add points from peakDoseForTinyMesh call
//     // pkSensiPlot->AddPoint(0.0001*1000., 4.13507);                  // doseForTinyMesh(0, 0.0001)
//     // pkSensiPlot->AddPoint(0.00015*1000., 1.88499);                 // doseForTinyMesh(0, 0.00015)
//     // pkSensiPlot->AddPoint(0.0002*1000., 1.12026);                  // doseForTinyMesh(0, 0.0002)
//     // pkSensiPlot->AddPoint(0.0003*1000., 0.539289);                 // doseForTinyMesh(0, 0.0003)
//     // pkSensiPlot->AddPoint(0.0004*1000., 0.297116);                 // doseForTinyMesh(0, 0.0004)
//     // pkSensiPlot->AddPoint(0.0005*1000., 0.21945);                  // doseForTinyMesh(0, 0.0005)
//     // pkSensiPlot->AddPoint(0.0006*1000., 0.152822);                 // doseForTinyMesh(0, 0.0006)
//     // pkSensiPlot->AddPoint(0.0007*1000., 0.110842);                 // doseForTinyMesh(0, 0.0007)
//     // pkSensiPlot->AddPoint(0.0008*1000., 0.0923902);                // doseForTinyMesh(0, 0.0008)
//     // pkSensiPlot->AddPoint(0.0009*1000., 0.0714925);                // doseForTinyMesh(0, 0.0009)
//     TCanvas* peakDoseSensitivityCanvas = new TCanvas("peakDoseSensitivity", "Sensitivity of the peak dose from the volume element");
//     peakDoseSensitivityCanvas->SetGridx(); peakDoseSensitivityCanvas->SetGridy();

//     int stat =1;
//     resetProgress();
//     for(float mSX=0.025; mSX < maxMesh; mSX+= (maxMesh-0.025)/points){
//         float status = (float)(stat)/(float)(points+1);
//         printProgress(status);

//         double pDY = peakDose(detNb, mSX);
//         // double pDY = doseForTinyMesh(detNb, mSX);
//         pkSensiPlot->AddPoint(mSX*1000., pDY*1000.);
//         stat++;
//     }
//     pkSensiPlot->Draw();


//     // TGraph* pkSensiPlotMethod2 = new TGraph();
//     // pkSensiPlotMethod2->SetLineColor(kRed);
//     // stat =1;
//     // resetProgress();
//     // for(float mSX=0.025; mSX < maxMesh; mSX+= (maxMesh-0.025)/points){
//     //     float status = (float)(stat)/(float)(points+1);
//     //     printProgress(status);

//     //     // double pDY = peakDose(detNb, mSX);
//     //     double pDY = doseForTinyMesh(detNb, mSX);
//     //     pkSensiPlotMethod2->AddPoint(mSX*1000., pDY*1000.);
//     //     stat++;
//     // }
//     // pkSensiPlotMethod2->Draw("same");

        
//     peakDoseSensitivityCanvas->Print(whichDetShort(detNb) + "_peakDSensMesh" + ".png");
//     return pkSensiPlot;
// }

// Compare the energy deposition spectrum dependence with respect to the production cut
// void compareEnergyCuts(int detNb=0){
//     loadAnotherFile("data/runXY_6_cut001_air10.root");
//     TH1F* cut001 = new TH1F(edepSpectrum(detNb));
//     totalEdep(detNb);

//     loadAnotherFile("data/runXY_6_cut01_air10.root");
//     TH1F* cut01 = new TH1F(edepSpectrum(detNb));
//     cut01->SetLineColor(kRed);
//     totalEdep(detNb);
 
//     loadAnotherFile("data/runXY_6_cut1_air10.root");
//     TH1F* cut1 = new TH1F(edepSpectrum(detNb));
//     cut1->SetLineColor(kYellow);
//     totalEdep(detNb);

//     loadAnotherFile("data/runXY_6_cut10_air10.root");
//     TH1F* cut10 = new TH1F(edepSpectrum(detNb));
//     cut1->SetLineColor(kGreen);
//     totalEdep(detNb);

//     TCanvas* compareCutsCanvas = new TCanvas("compareCuts", "Edep comparison with cuts");
//     cut001->Draw();
//     cut01->Draw("same");
//     cut1->Draw("same");
//     cut10->Draw("same");
//     compareCutsCanvas->Print(whichDetShort(detNb) + "_edepCutDep" + ".png");
// }

// Compare the total/peak dose dependence with respect to the production cut
// void compareDoseCut(int detNb=0){
//     TGraph* doseComparePlot = new TGraph();
//     doseComparePlot->SetTitle("Dose sensitivity with respect to cut value");
//     doseComparePlot->GetXaxis()->SetTitle("Cut value [um]");
//     doseComparePlot->GetYaxis()->SetTitle("Total dose [mGy]");
    
//     TGraph* peakdoseComparePlot = new TGraph();
//     peakdoseComparePlot->SetTitle("Peak dose sensitivity with respect to cut value");
//     peakdoseComparePlot->GetXaxis()->SetTitle("Cut value [um]");
//     peakdoseComparePlot->GetYaxis()->SetTitle("Peak dose [mGy]");

//     loadAnotherFile("data/runXY_7_cut001_air10.root");
//     doseXY(detNb);
//     doseComparePlot->AddPoint(0.01, totDoseGlobal*1E3);
//     peakdoseComparePlot->AddPoint(0.01, peakDoseGlobal*1E3);


//     loadAnotherFile("data/runXY_7_cut01_air10.root");
//     doseXY(detNb);
//     doseComparePlot->AddPoint(0.1, totDoseGlobal*1E3);
//     peakdoseComparePlot->AddPoint(0.1, peakDoseGlobal*1E3);
//     // cut01->SetLineColor(kRed);
 
//     loadAnotherFile("data/runXY_7_cut1_air10.root");
//     doseXY(detNb);
//     doseComparePlot->AddPoint(1, totDoseGlobal*1E3);
//     peakdoseComparePlot->AddPoint(1, peakDoseGlobal*1E3);
//     // cut1->SetLineColor(kYellow);

//     loadAnotherFile("data/runXY_7_cut10_air10.root");
//     doseXY(detNb);
//     doseComparePlot->AddPoint(10, totDoseGlobal*1E3);
//     peakdoseComparePlot->AddPoint(10, peakDoseGlobal*1E3);
//     // cut1->SetLineColor(kGreen);

//     loadAnotherFile("data/runXY_7_cut100_air10.root");
//     doseXY(detNb);
//     doseComparePlot->AddPoint(100, totDoseGlobal*1E3);
//     peakdoseComparePlot->AddPoint(100, peakDoseGlobal*1E3);

//     loadAnotherFile("data/runXY_7_cut1000_air10.root");
//     doseXY(detNb);
//     doseComparePlot->AddPoint(1000, totDoseGlobal*1E3);
//     peakdoseComparePlot->AddPoint(1000, peakDoseGlobal*1E3);

//     TCanvas* doseSensitivityCanvas = new TCanvas("doseSensitivityCanvas", "Sensitivity of the total dose from the volume element");
//     doseSensitivityCanvas->SetGridx(); doseSensitivityCanvas->SetGridy();
//     doseComparePlot->SetMarkerStyle(30);
//     doseComparePlot->GetXaxis()->SetLimits(0.01-0.2,1000.0+1.);
//     doseSensitivityCanvas->SetLogx();
//     doseComparePlot->Draw();
//     doseSensitivityCanvas->Print(whichDetShort(detNb) + "_doseSensCut.png");

//     TCanvas* peakdoseSensitivityCanvas = new TCanvas("peakdoseSensitivityCanvas", "Sensitivity of the peak dose from the volume element");
//     peakdoseSensitivityCanvas->SetGridx(); peakdoseSensitivityCanvas->SetGridy();
//     peakdoseComparePlot->SetMarkerStyle(30);
//     peakdoseComparePlot->GetXaxis()->SetLimits(0.01-0.2,1000.0+1.);
//     peakdoseSensitivityCanvas->SetLogx();
//     peakdoseComparePlot->Draw();
//     peakdoseSensitivityCanvas->Print(whichDetShort(detNb) + "_dosePkSensCut.png");
// }

// // compare ?? duplicate?
// void compareDoseCuts(int detNb=0){
//     float cutList[] = {10., 1., 0.1, 0.01};
//     TString cutListText[] = {"10", "1", "01", "001"};
//     float yPoints[4];

//     TGraph* dsComparePlot = new TGraph();
//     dsComparePlot->SetTitle("Dose sensitivity with cut value");
//     dsComparePlot->GetXaxis()->SetTitle("Cut value [um]");
//     dsComparePlot->GetYaxis()->SetTitle("Dose [Gy]");

//     // std::to_string(cutList[i])
    
//     for(int i=0; i< size(cutList); i++){
//         TString filename = "data/runXY_cut" + cutListText[i] + ".root";
//         cout << filename << endl;
//         loadAnotherFile(filename);
//         doseXY(detNb);
//         dsComparePlot->AddPoint(cutList[i], peakDoseGlobal);
//         yPoints[i] = peakDoseGlobal;
//     }

//     // size(cutList), cutList, yPoints

//     TCanvas* dsComparePlotCanvas = new TCanvas("dsComparePlotCanvas", "Sensitivity of the total dose from the production cut");
//     dsComparePlotCanvas->SetGridx(); dsComparePlotCanvas->SetGridy();
//     dsComparePlot->GetXaxis()->SetLimits(0.1-0.5,10.5);
//     dsComparePlot->SetMarkerStyle(30);
    
//     dsComparePlot->Draw();
//     dsComparePlotCanvas->Print(whichDetShort(detNb) + "_peakDSens.png");
// }