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


// TString ifilename = "kyle/kyleBeam_6_50_0_w100.root";
TString ifilename = "MIP.root";

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
    std::vector<TH1D> beamMonitor(int bins = 200, float range = 10);                            // Primary particle monitor
    std::vector<TH1D> beamMonitor(int pdg, int bins = 200, float range = 10);                   // Primary particle monitor for initial particles with assigned pdg code
    std::vector<TH2D> primaryProfile(int pdg=11);                                               // Number XY profile of the particles. Filter for particles with given pdg code
    std::vector<TH1D> primarySpectrum();                                                        // Energy spectrum of primaries
    std::vector<TH1D> primarySpectrum(int pdg, double range = 10.0);                            // Energy spectrum of primaries with the given pdg code. Return 2-array with 0-pdgfiltered, 1-all part.s
    std::vector<TH1D> primarySpectrumLogLog(double lowRange = 1E-6,                             // Energy spectrum of primaries with the given pdg code, with variable size bins. Return 2-array with 0-pdgfiltered, 1-all part.s
        double highRange = 16.0, int binsNb = 100);
    std::vector<TH2D> primaryEnergyProfile();                                                   // Energy XY profile of the particles. Each point (x,y) is weightned with the corresponding energy
// Energy deposition plots
    std::vector<TH1D> edepSpectrum(int bins, float emin, float emax);                           // Energy spectrum of depositions for the upstream/downstream detector
    std::vector<TH1D> edepSpectrumPartitionedByPrimaryEnergies(int bins=1000,                   // Energy spectrum of depositions for the upstream/downstream detector in a window Ea-Eb, partitioned by primary energies. Energy in keV
    double emin = 0.0 /*GeV*/, double emax= 6.0 /*GeV*/, double baseRatio = -1E6 /*keV*/);
    std::vector<TH1D> edepSpectrumChg(int bins, double a, double b);                            // Impact of charged particles to the spectrum of deposited energy. Useful to understand contribution to the spectrum from gamma conversions in different areas of the detector. Also for MIP. Return a 4-array with the spectrum of energy deposited in upstream/downstream detector in the range [a,b] with #bins. Sel histograms are filled only with events where a incoming pdg !=22 track enters. [0] det0 [1] det1 [2] det0_sel [3] det1_sel
    std::vector<TH1D> edepSpectrumChg();                                                        // Overloading
    std::vector<TH1D> eSpectrumStrip(int stripNb = 100, float range = 5.);                      // Energy spectrum of depositions in a given strip number for some det#
    double totalEdep(int detNb=0);                                                              // Energy deposited in detNb in keV
    double* eStrip_array(int stripNb = 100);                                                    // 2-array of energy deposited in the strip 100 for upstream/downstream detectors in keV
    double eStrip(int detNb=0, int stripNb = 100);                                              // Energy deposited in the strip 'stripNb' of detector 'detNb' in keV
    // double eStrip(int detNb=0, int stripNb = 100, float range = 1.);                         // Energy spectrum of depositions in a given strip number. Return energy deposited in keV
    std::vector<TH1D> eLongitudinal_array(int nbBins = 100);                                    // Longitudinal energy deposition profile. Custom binning available. Return 2-array with det0 and det1 plots
    TH1D eLongitudinal(int detNb);                                                              // Longitudinal energy deposition profile for detNb
    TH1D eLongitudinal(int detNb, int nbBins);                                                  // Longitudinal energy deposition profile for detNb with custom binning.
    std::vector<TH1D> edepProject(double meshSize = 0.100);                                     // Energy deposition in the transverse plane. Return 4-array [0]det0x, [1]det0y, [2]det1x, [3]ddet1y
    std::vector<TH1D> edepProject(bool filter, double meshSize = 0.100, double norm=1);         // Energy deposition in the transverse plane for primaries with a selection filter (i.e. charge/neutral). Return 4-array [0]det0x, [1]det0y, [2]det1x, [3]ddet1y
    std::vector<TH1D> edepProjectStripTree(double norm=1);                                      // Energy deposition profile in the strips. It uses the strip tree data. Return a 2-array for det0 and 1
    std::vector<TH1D> edepProjectStripTree(double pE0, double pE1, double norm=1);              // Energy deposition profile in the strips given primaries in selected energy range [pE0, pE1] GeV. Return a 2-array for det0 and 1
    std::vector<TH1D> edepStripChg(double norm = 1);                                            // Energy deposition profile in the strips with comparison with deposition from charged primaries, for the primaries in the energy range [enLow, enHigh]. It uses the strip tree data. Return a 2-array for det0 and 1
    std::vector<TH1D> edepProjectPrimEnPartitioned(double enLow=0 /*GeV*/,                      // Energy depositions in the strips separated between the contributions with different energies
        double enHigh=10.0 /*GeV*/, double baseRatio = 10.0);
    void plot_edepProjectPrimEnPartitioned(double enLow=0 /*GeV*/,                              //TEMPDESCRIPTION-Plot function edepProjectPrimEnPartitioned
        double enHigh=10.0 /*GeV*/, double baseRatio = 10.0, 
        TString path = "plots/",TString format = ".pdf");
// Dose
    std::vector<TH2D> doseXY(double meshSize = 0.100 /*mm*/);                                   // Energy & dose 2D map using a mesh of size LxLx100 um3. Return a 4-array of TH2D with [0]emapdet0, [1]emapdet1, [2]dmapdet0, [3]dmapdet1
    std::vector<TH2D> doseXY(int filter, double meshSize = 0.100 /*mm*/);                       // Energy & dose 2D map using a mesh of size LxLx100 um3. Contributions from only photons to the dose. Return a 4-array of TH2D with [0]emapdet0, [1]emapdet1, [2]dmapdet0, [3]dmapdet1
// Profile reconstruction    
    std::vector<TH1I> fitHits(double threshold, int* fitRange);                                 // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Threshold is in keV. Return a 2-array with fit for upstream and downstream det.s
    std::vector<TH1I> fitHits();                                                                // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Overloading
    double* fitHitsChi2(double threshold, int* fitRange);                                       // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Threshold is in keV, fitRange is a 2-array with stripL, stripH where the fit is ranged.
    double* fitHitsChi2(double threshold = 1.);                                                 // Energy depositions/det over X (Y) & profile reconstruction via the 'threshold method'. Default threshold of 1keV and range 90,110.
    std::vector<TGraph2D> fitVSThrRng();                                                        // 2D-graph of the chisquare as a function of threshold and range
    void plot_fitVSThrRng();                                                                    //TEMPDESCRIPTION-Plot fitVSThrRng()
// Other methods + unsorted
    void overlayReconstructedBeamProfile();                                                     // Superimpose the beam profile of the primaries with the reconstructed one
    // int createReport(TString path = "plots/", TString format = ".pdf",                          // Save plots in the folder plots
    //     TString ofilename = "out.root");
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

// Rescale to the BX
double rescaleBX(){
    int config = getConfig(); double factor = 0;
    if(config > 0 || (config == -1 && ifilename.Contains("kyleBeam"))){
        factor = pow(10.0, 9.0 - config);
    }else if(config == -1){
        cout << "Config not found. getConfig() = " << config << endl;
    }else{
        factor = -1500.0/config;
    }
    return factor;
}



// Primary particle monitor
std::vector<TH1D> beamMonitor(int bins = 200, float range = 10){
    TH1D result[3];
    result[0] = TH1D("beamX", "Beam monitor;x [mm];counts", bins, -range, range);
    result[1] = TH1D("beamY", "Beam monitor;y [mm];counts", bins, -range, range);
    result[2] = TH1D("beamZ", "Beam monitor;z [mm];counts", 110, -110, 0);
    primaryTree->Project("beamX", "x0");
    primaryTree->Project("beamY", "y0");
    primaryTree->Project("beamZ", "z0");

    // Rescale to the BX
    double factor = rescaleBX();
    result[0].Scale(factor);
    result[1].Scale(factor);
    result[2].Scale(factor);

    if(verbosityLevel > 0){
        TCanvas* beamMonitorCanvas = new TCanvas("beamMonitorCanvas", "Primary particle monitor", 1600, 400);
        beamMonitorCanvas->Divide(3);
        for(int i=0; i<3; i++){
            beamMonitorCanvas->cd(i+1);
            result[i].DrawClone("hist");
        }
    }
    std::vector<TH1D> retVar;
    for(int i=0; i<3; i++){
        retVar.push_back(result[i]);
    }
    return retVar;
}

// Primary particle monitor for initial particles with assigned pdg code
std::vector<TH1D> beamMonitor(int pdg, int bins = 200, float range = 10){
    TH1D result[3];
    result[0] = TH1D("beamX", "Beam monitor;x [mm];counts", bins, -range, range);
    result[1] = TH1D("beamY", "Beam monitor;y [mm];counts", bins, -range, range);
    result[2] = TH1D("beamZ", "Beam monitor;z [mm];counts", 110, -110, 0);
    primaryTree->Project("beamX", "x0", TString::Format("pdg==%i", pdg));
    primaryTree->Project("beamY", "y0", TString::Format("pdg==%i", pdg));
    primaryTree->Project("beamZ", "z0", TString::Format("pdg==%i", pdg));

    // Rescale to the BX
    double factor = rescaleBX();
    result[0].Scale(factor);
    result[1].Scale(factor);
    result[2].Scale(factor);

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

    // Rescale to BX
    double factor = rescaleBX();
    primaryeSpectrHist->Scale(factor);
    primaryeSpectrPDGHist->Scale(factor);

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
std::vector<TH1D> primarySpectrumLogLog(double lowRange = 1E-6, double highRange = 16.0, int binsNb = 100){
    // Creating the bins array 
    double* bins = (double*)malloc((binsNb+1)*sizeof(double));
    bins[0] = lowRange;
    double N = (binsNb)/( 1 - (TMath::Log(lowRange)/TMath::Log(highRange)));
    for(int i=0; i<binsNb; i++){
        bins[binsNb-i] = pow(highRange, 1 - i/N);
        // cout << bins[binsNb-i] << '\t';
    }
    // cout << bins[0] << endl;

    // cout << "breakpoint A. pars: " << pdg << '\t' << range << endl;
    TH1D* primaryeSpectrHist = new TH1D("primaryeSpectrLogBinsHist", "Energy spectrum of primaries; energy [GeV]; counts", binsNb, bins);
    TH1D* primaryeSpectrPDGHist = new TH1D("primaryeSpectrPDGLogBinsHist", "Energy spectrum of primaries with pdg != 22; energy [GeV]; counts", binsNb, bins);
    // cout << "breakpoint A2. pars: " <<  range*10.0 << '\t' <<  0 << '\t' << range << endl;
    primaryTree->Project("primaryeSpectrLogBinsHist", "ekin");
    // cout << "breakpoint A22. pars: " <<  range*10.0 << '\t' <<  0 << '\t' << range << endl;
    primaryTree->Project("primaryeSpectrPDGLogBinsHist", "ekin", "pdg!=22");
    // cout << "breakpoint A3" << endl;
    primaryeSpectrPDGHist->SetLineColor(kRed);
    primaryeSpectrHist->GetXaxis()->SetRangeUser(lowRange, 1.5*highRange);
    primaryeSpectrPDGHist->GetXaxis()->SetRangeUser(lowRange, 1.5*highRange);

    // Rescale to BX
    double factor = rescaleBX();
    primaryeSpectrHist->Scale(factor);
    primaryeSpectrPDGHist->Scale(factor);
    primaryeSpectrHist->Scale(factor);
    primaryeSpectrPDGHist->Scale(factor);


    if(verbosityLevel>0){
        TCanvas* primarySpectrumPDGCanvas = new TCanvas("primarySpectrumPDGCanvas", "Energy spectrum of pdg!=22 primaries");
        primarySpectrumPDGCanvas->SetGridx(); primarySpectrumPDGCanvas->SetGridy();
        primarySpectrumPDGCanvas->SetLogx(); primarySpectrumPDGCanvas->SetLogy();
        TLegend* legend = new TLegend(0.25,0.8,0.45,0.9);
        legend->SetBorderSize(0);   legend->SetFillColor(0);    legend->SetTextSize(0.04);
        legend->AddEntry(primaryeSpectrHist, "#gamma, others", "l");
        legend->AddEntry(primaryeSpectrPDGHist, "others", "l");
        primaryeSpectrHist->DrawClone("hist");
        primaryeSpectrPDGHist->DrawClone("hist same");
        legend->Draw();
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
// Impact of charged particles to the spectrum of deposited energy. Useful to understand contribution to the spectrum from gamma conversions in different areas of the detector. Also for MIP. Return a 4-array with the spectrum of energy deposited in upstream/downstream detector in the range [a,b] with #bins. Sel histograms are filled only with events where a incoming pdg !=22 track enters. [0] det0 [1] det1 [2] det0_sel [3] det1_sel
std::vector<TH1D> edepSpectrumChg(int bins = 500, double a=0, double b=500){
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
        if(pdgSelU != 22){
            pSelListU.push_back(eventSelU);
        }
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the upstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    // std::sort(pSelListU.begin(), pSelListU.end());

    
    
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
    std::sort(pSelListD.begin(), pSelListD.end(),
        // Lambda expression begins
        [](pSelListDCl a, pSelListDCl b) {
            return (a.getTrack() > b.getTrack());
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction

    // Dose code with some minor modifications
    double totEnergy[2] = {0};

    TH1D* eSpectrHist[2];
    TH1D* eSpectrSelHist[2];
    for(int i=0; i<2; i++){
        eSpectrHist[i] = new TH1D(TString::Format("edepSpectrDet%i",i), "Energy spectrum " + whichDetVerbose(i) + " detector;keV; counts", bins, a, b);
        eSpectrSelHist[i] = new TH1D(TString::Format("edepSpectrSelDet%i",i), "Energy spectrum " + whichDetVerbose(i) + " detector pdg!=22;keV; counts", bins, a, b);
        // Style
        eSpectrSelHist[i]->SetLineColor(kRed);
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

    int prevEvent=1;
    double edepInEvent[2] = {};
    double edepInEventSel[2] = {};

    int pU0=0; int pD0=0;
    resetProgress();
    for(Long64_t i=0; i< debugEntries; i++ ){
        debugTree->GetEntry(i);
        float status = (float)(i+1)/(float)debugEntries;
        printProgress(status);

        if(det >= 2) continue; // Ignore the charge scorer plane
        totEnergy[det] += edep;
        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\t";

        // The event is the same of the previous one. Keep summing energy dep. contribution in the same arrays (edepInEvent or edepInEventSel)
        if(event == prevEvent){
            // cout << "=";
            // Energy contributions from all tracks without any tagging
            edepInEvent[det] += edep;

            // Energy contribution from selected tracks
            if(det==0){
                // cout << " det0\n";
                // For det0 selection is made using eventID since only 1 track is injected in sapphire because of the way the beam is imported
                for(int j=pU0; j<selUNb; j++){
                    int val=pSelListU[j];
                    if(val > event){
                        break;
                    }else if(val == event){
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\tchg\n";
                        edepInEventSel[det] += edep;
                        pU0=j;
                        break;
                    }
                }
            }else if(det==1){ // <- I know, but this way it's clearer...
                // cout << " det1\n";
                // For det1 selection is more tricky. It is made using track IDs of the scoring plane. So a for loop must be done.
                for(int j_events=pD0; j_events<selDNb; j_events++){
                    auto pSelClass_ith = pSelListD[j_events];
                    int evtOf_ith_class = pSelClass_ith.getEvent();
                    int trkOf_ith_class = pSelClass_ith.getTrack();
                    if(event < evtOf_ith_class){
                        break;
                    }else if(event == evtOf_ith_class){
                        if(track == trkOf_ith_class){
                            // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\tchg\n";
                            edepInEventSel[det] += edep;
                            pD0=j_events;
                            break;
                        }
                        continue;
                    }
                }
            }
        }else{
            // cout << "- ";
            // The event changed. Store the value for the edep and continue
            for(int ii=0; ii<2; ii++){
                // cout << edepInEvent[ii] << ", " << edepInEventSel[ii] << "; ";
                if(edepInEvent[ii] != 0){
                    eSpectrHist[ii]->Fill(edepInEvent[ii]);
                    // cout << "filling htot with: " << edepInEvent[ii] << endl;
                    edepInEvent[ii] = 0;
                }
                if(edepInEventSel[ii] != 0){
                    eSpectrSelHist[ii]->Fill(edepInEventSel[ii]);
                    // cout << "filling hSel with: " << edepInEventSel[ii] << endl;
                    edepInEventSel[ii] = 0;
                }
            }

            // Add this entry deposition
            edepInEvent[det] += edep;
            // Energy contribution from selected tracks
            if(det==0){
                // cout << " det0\n";
                // For det0 selection is made using eventID since only 1 track is injected in sapphire because of the way the beam is imported
                for(int j=pU0; j<selUNb; j++){
                    int val=pSelListU[j];
                    if(val > event){
                        break;
                    }else if(val == event){
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\tchg\n";
                        edepInEventSel[det] += edep;
                        pU0=j;
                        break;
                    }
                }
            }else if(det==1){
                // cout << " det1\n";
                // For det1 selection is more tricky. It is made using track IDs of the scoring plane. So a for loop must be done.
                for(int j_events=pD0; j_events<selDNb; j_events++){
                    auto pSelClass_ith = pSelListD[j_events];
                    int evtOf_ith_class = pSelClass_ith.getEvent();
                    int trkOf_ith_class = pSelClass_ith.getTrack();
                    if(event < evtOf_ith_class){
                        break;
                    }else if(event == evtOf_ith_class){
                        if(track == trkOf_ith_class){
                            // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\tchg\n";
                            edepInEventSel[det] += edep;
                            pD0=j_events;
                            break;
                        }
                        continue;
                    }
                }
            }
            prevEvent = event;
        }
    }

    if(verbosityLevel>0){
        cout << "tot edep in det0: " << totEnergy[1] << endl; 
        TCanvas* edepSpectrumCanvas = new TCanvas("edepSpectrumCanvas", "Energy spectrum", 1200, 600);
        edepSpectrumCanvas->Divide(2,2);
        edepSpectrumCanvas->cd(1);
        eSpectrHist[0]->Draw();
        edepSpectrumCanvas->cd(2);
        eSpectrHist[1]->Draw();
        edepSpectrumCanvas->cd(3);
        eSpectrSelHist[0]->Draw();
        edepSpectrumCanvas->cd(4);
        eSpectrSelHist[1]->Draw();
    }

    std::vector<TH1D> retVar = {TH1D(*eSpectrHist[0]), TH1D(*eSpectrHist[1]), TH1D(*eSpectrSelHist[0]), TH1D(*eSpectrSelHist[1]) };
    return retVar;
}
// Overload
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
std::vector<TH1D> edepProjectStripTree(double norm = 1){
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
    
    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        result[i]->Scale(factor * 1E-6);
        result[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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
std::vector<TH1D> edepProjectStripTree(double pE0, double pE1, double norm = 1){
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

    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        edepStripSelHist[i]->Scale(factor * 1E-6);
        edepStripSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        edepStripHist[i]->Scale(factor * 1E-6);
        edepStripHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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
std::vector<TH1D> edepStripChg(double norm = 1){
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1;
    std::vector<TH1D> tmp = edepProject(0, 0.1, norm);
    verbosityLevel = verbosityLevel_bak;
    std::vector<TH1D> retVar = {tmp[0], tmp[3], tmp[4], tmp[7]};
    for(int i=0; i<4; i++){
        retVar[i].GetXaxis()->SetLimits(1,200);
        retVar[i].GetXaxis()->SetTitle("strip number");
    }
    retVar[0].SetTitle("Energy deposited per strip in the upstream detector");
    retVar[1].SetTitle("Energy deposited per strip in the downstream detector");
    retVar[2].SetTitle("Energy deposited per strip in the upstream detector");
    retVar[3].SetTitle("Energy deposited per strip in the downstream detector");

    if(verbosityLevel>0){
        TCanvas* edepStripChgCanvas = new TCanvas("edepStripChgCanvas", "edepStripChgCanvas", 1000, 400);
        edepStripChgCanvas->Divide(2);
        edepStripChgCanvas->cd(1);
        retVar[0].DrawClone("hist");
        retVar[2].DrawClone("hist same");
        edepStripChgCanvas->cd(2);
        retVar[1].DrawClone("hist");
        retVar[3].DrawClone("hist same");
    }
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
    const int nbYpts = 20./meshSize;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector BX [GeV];X [mm]; Y [mm]", nbXpts, -10.0, 10.0, nbYpts, -10.0, 10.0);
        // edepXYHist[i].GetXaxis()->SetRangeUser(-1.5, 1.5);
        // edepXYHist[i].GetYaxis()->SetRangeUser(-1.5, 1.5);
        // Style
        // edepXYHist[i]->SetStats(0);
        // edepXYHist[i]->SetContour(500);
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
        if(det == 2) continue;
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

    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        edepXYHist[i]->Scale(factor * 1E-6);
        doseXYHist[i]->Scale(factor);
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
            if(verbosityLevel>1){
                cout    << "Total energy deposit is:\t"             << totEnergy[i] * 1E-6 << " GeV" << endl;
                cout    << "Total absorbed dose is:\t\t"            << totEnergy[i]/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
                cout    << "Peak dose is:\t\t\t"                    << peakDose[i] << " Gy (" << peakEnergy[i]*1E3 << " MeV in " << vol << " um3)" <<  endl;
                cout    << "----------------BX--------------------" << endl;
            }
            cout    << "Total energy deposit/BX is:\t"          << totEnergy[i] * factor * 1E-6 << " GeV" << endl;
            cout    << "Total absorbed dose/BX is:\t"           << (totEnergy[i] * factor)/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
            cout    << "Peak dose/BX is:\t\t"                   << peakDose[i] << " Gy (" << peakEnergy[i] * 1E3 << " MeV in " << vol << " um3)" <<  endl;
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
std::vector<TH2D> doseXY(int filter, double meshSize = 0.100 /*mm*/){
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
        if(pdgSel != filter) pSelList.push_back(i);
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
        // edepXYHist[i]->SetContour(500);
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

    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        edepXYHist[i]->Scale(factor * 1E-6);
        doseXYHist[i]->Scale(factor);
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
            if(verbosityLevel>1){
                cout    << "Total energy deposit is:\t"             << totEnergy[i] * 1E-6 << " GeV" << endl;
                cout    << "Total absorbed dose is:\t\t"            << totEnergy[i]/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
                cout    << "Peak dose is:\t\t\t"                    << peakDose[i] << " Gy (" << peakEnergy[i]*1E3 << " MeV in " << vol << " um3)" <<  endl;
                cout    << "----------------BX--------------------" << endl;
            }
            cout    << "Total energy deposit/BX is:\t"          << totEnergy[i] * factor * 1E-6 << " GeV" << endl;
            cout    << "Total absorbed dose/BX is:\t"           << (totEnergy[i] * factor)/totVol * conversionFactor << " Gy (" << totVol << " um3)" << endl;
            cout    << "Peak dose/BX is:\t\t"                   << peakDose[i] << " Gy (" << peakEnergy[i] * 1E3 << " MeV in " << vol << " um3)" <<  endl;
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
// Energy deposition in the transverse XY plane. Either X or Y profile by selecting the returned value [0] det0X [1] det0Y [2] det1X [3] det1Y
std::vector<TH1D> edepProject(bool filter, double meshSize = 0.100, double norm=1){
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
        if(pdgSelU != 22){
            pSelListU.push_back(eventSelU);
        }
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the upstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListU.begin(), pSelListU.end());

    
    
    // Downstream detector.
    // Selection happens using the track ID and the charge scoring plane tree
    class pSelListDCl {
        public:
        pSelListDCl(int event, int track){ eventID = event; trackID = track;};
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
            // cout << trackSelD << ", ";
            pSelListDCl tmp(eventSelD, trackSelD);
            pSelListD.push_back(tmp);
        }
    }
    Long64_t selDNb = pSelListD.size();
    if(verbosityLevel>0) cout << "Non-photons entering the downstream detector: " << selDNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListD.begin(), pSelListD.end(),
        // Lambda expression begins
        [](pSelListDCl a, pSelListDCl b) {
            return (a.getEvent() < b.getEvent() && a.getTrack() < b.getTrack());
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction
    // Debug
    // for(int i=0; i<selDNb; i++){
    //     cout << pSelListD[i].getEvent() << ":" <<  pSelListD[i].getTrack() << "\t";
    // }
    
    // Dose code with some minor modifications
    const int binsNb = 20./meshSize;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    TH2D* edepXYSelHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        edepXYSelHist[i] = new TH2D(TString::Format("edepXYSelHist%i",i), "non-photons primaries " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
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
                }else if(val == event){
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
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\t-chg\n";
                        edepXYSelHist[det]->Fill(x, y, edep);
                        pD0=j_events;
                        break;
                    }
                    continue;
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


    // Rescale to BX
    double factor = rescaleBX();
    if(factor == 0 && norm != 1){
        cout << "Manual normalization detected. This file corresponds to " << norm << " BX. Rescaling accordingly..." << endl;
        factor = norm * 1E-6;
    }
    
    //
    for(int i=0; i<2; i++){
        edepXYHist[i]->Scale(factor * 1E-6);
        edepXYSelHist[i]->Scale(factor * 1E-6);
        edepXYHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        edepXYSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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
        edepProj[i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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































// Energy deposition in the transverse XY plane. Either X or Y profile by selecting the returned value [0] det0X [1] det0Y [2] det1X [3] det1Y
std::vector<TH1D> edepProjectMOOOOD(bool filter, double meshSize = 0.100, double norm=1){
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
        if(pdgSelU != 22){
        // if(1){
            pSelListU.push_back(eventSelU);
        }
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the upstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListU.begin(), pSelListU.end());

    
    // for(int i = 0; i<pSelListU.size(); i++){
    //     cout << pSelListU[i] << ", ";
    // }
    // cout << endl << "------------------\n\n\n";


    // Pop events that correspond to particles created in the air between the kapton and the GBP
    std::vector<int> exclusionTable = {1038, 1039, 1667, 1668, 4507, 4508, 5397, 5398, 5440, 5913, 5914, 6319, 6320, 6780, 8562, 8563, 11637, 13887, 13888, 15801, 18282, 18283, 18412, 18413, 19917, 19918, 20514, 20515, 27212, 27213, 27439, 27440, 30168, 30169, 30575, 31376, 31377, 32064, 32065, 34133, 34134, 35729, 35730, 35878, 35879, 38687, 38688, 39372, 39373, 39981, 39982, 41138, 42606, 42607, 43202, 43203, 46990, 48513, 48514, 48721, 48722, 48937, 49163, 49164, 49923, 49924, 51799, 51800, 52495, 53021, 53022, 54911, 54912, 55281, 55921, 55922, 56282, 57813, 57814, 57843, 57844, 58496, 58497, 58541, 58542, 58793, 58794, 58964, 58965, 61279, 61280, 61410, 61411, 62974, 63054, 63055, 64433, 64434, 66874, 66875, 70445, 70446, 70552, 70553, 71032, 71033, 73332, 73333, 73930, 74478, 74479, 74823, 74824, 75857, 77148, 77149, 78940, 78941, 79644, 79645, 80603, 81904, 83447, 83448, 85120, 85121, 86580, 86581, 87081, 88169, 88170, 92390, 92391, 94032, 94033, 98726, 98727, 101753, 101754, 106061, 106062, 108916, 108917, 110209, 110210, 111058, 111059, 114285, 114286, 115854, 115855, 118372, 118373, 119380, 119381, 120600, 120601, 121487, 121488, 122014, 122015, 127153, 127283, 127284, 128010, 128011, 129361, 129362, 130100, 130101, 130513, 130514, 133155, 133156, 135648, 135649, 138986, 138987, 139295, 140152, 141432, 141433, 146493, 151703, 151704, 152330, 153291, 153292, 153873, 153874, 159054, 159055, 160805, 160806, 162166, 162167, 162272, 162273, 164466, 164467, 166464, 166465, 167230, 167338, 167339, 168038, 168039, 173015, 173016, 175550, 175551, 178753, 178754, 183003, 183004, 183181, 183858, 184492, 184493, 184794, 184795, 185619, 185620, 185657, 186005, 186006, 186198, 186199, 188868, 188869, 189270, 189271, 189535, 189536, 190938, 192664, 192665, 194383, 194384, 195255, 195256, 198078, 198079, 200026, 200027, 201298, 201299, 202674, 202675, 203935, 203936, 206662, 206663, 207327, 207328, 209823, 209846, 209847, 211781, 211782, 211897, 212035, 212036, 212324, 212325, 214905, 214906, 215887, 217093, 217094, 218312, 218313, 218432, 220831, 220832, 221950, 221951, 222568, 222569, 222897, 222898, 224812, 224813, 225028, 225029, 225556, 225557, 226267, 226268, 227267, 227268, 227949, 227950, 232194, 232195, 233754, 233755, 235020, 235021, 235181, 235182, 235490, 237206, 237234, 237235, 237343, 238407, 239513, 239514, 241353, 241354, 242337, 242338, 244542, 244543, 244980, 244981, 248730, 248731, 249013, 249014, 251546, 251547, 252324, 252325, 253347, 253348, 253936, 253937, 255330, 255331, 255655, 260586, 260587, 261230, 261231, 264834, 264835, 266867, 266868, 267704, 267705, 269605, 270798, 270799, 270827, 270828, 271496, 271497, 273542, 273543, 276252, 276253, 276993,
        276994, 278162, 278163, 283267, 283268, 283779, 283780, 284210, 284211, 285497, 285498, 285887, 287332, 287333, 289239, 289321, 289322, 296867, 296913, 296914, 297649, 297650, 298033, 298034, 299662, 299663, 300408, 300409, 303123, 303124, 303558, 303559, 306134, 306135, 306924, 306925, 307311, 307312, 307555, 307556, 309673, 310918, 310919, 311068, 311069, 311995, 311996, 314647, 314648, 314686, 314687, 315364, 315365, 315577, 315578, 315695, 315696, 316933, 316934, 317221, 317222, 318101, 318102, 318103, 318123, 318124, 318631, 318632, 319764, 319765, 320833, 320834, 321342, 321343, 321884, 321885, 325091, 325092, 326167, 326168, 327158, 327159, 327236, 327237, 327443, 327444, 330063, 330064, 330939, 330940, 330941, 333360, 333361, 333688, 333689, 334736, 334737, 334987, 334988, 335219, 335220, 338176, 338177, 338178, 338179, 339395, 339396, 340486, 340487, 342003, 342004, 344288, 344457, 344458, 344603, 344604, 347610, 347611, 347955, 347956, 347967, 347968, 349181, 349182, 350548, 350549, 350620, 350621, 350925, 350926, 352057, 352058, 352502, 352503, 352665, 352666, 353111, 353112, 355570, 355571, 356445, 356446, 356562, 356973, 356974, 357220, 357221, 358103, 358104, 358652, 358653, 358985, 358986, 361292, 361293, 363375, 363376, 363821, 365084, 365085, 365987, 367810, 367978, 371341, 371342, 375789, 375790, 377320, 377321, 377651, 377652, 377763, 377764, 382099, 382100, 382739, 382740, 382777, 382778, 383807, 383808, 384604, 385896, 385897, 385904, 385905, 386035, 394647, 394648, 396555, 396556, 397331, 397332, 398297, 398395, 398396, 398630, 398631, 400249, 400250, 402147, 402148, 402449, 402450, 402534, 402554, 402555, 409291, 409292, 409571, 409572, 411074, 411075, 412806, 412807, 416216, 416217, 416988, 416989, 418331, 418332, 419377, 419378, 423265, 423266, 424685, 424686, 426621, 426622, 426834, 426835, 427738, 427739, 430815, 430816, 432125, 432126, 433121, 433122, 433401, 433402, 434000, 434001, 438100, 438729, 438730, 438942, 438943, 439982, 439983, 440666, 440668, 442589, 442590, 443496, 443497, 444042, 444043, 444121, 444122, 447677, 447678, 448660, 448661, 449829, 449830, 453990, 453991, 454062, 454107, 455619, 461611, 461612, 465074, 465075, 465886, 465887, 466130, 466131, 466132, 467154, 467155, 469829, 469830, 472213, 472214, 472283, 472284, 475226, 475227, 475835, 475836, 476939, 476940, 477737, 477738, 478832, 478833, 479542, 479543, 482036, 482037, 483594, 483595, 484485, 484486, 485185, 485186, 486606, 486607, 486886, 489597, 489598, 494732, 494733, 495080, 495081, 495179, 496326, 496327, 501144, 501678, 501679, 503348, 503349, 504967, 504968, 505265, 505266, 505873, 506016, 506017, 508393, 508394, 508582, 508583, 508848, 508849, 512692,
        512693, 513491, 513492, 514198, 514199, 517747, 517748, 517913, 517914, 517947, 517948, 518135, 518136, 519654, 520489, 520490, 523112, 523113, 524359, 524360, 525711, 525712, 525848, 525849, 529759, 529760, 530852, 530853, 533298, 533299, 533726, 533727, 536083, 536084, 540271, 540272, 540832, 540833, 540989, 540990, 542955, 542956, 543950, 543951, 545745, 545746, 546048, 546730, 559883, 559884, 560383, 560384, 562568, 566109, 566110, 566677, 566678, 567982, 567983, 569416, 569417, 569667, 569668, 571867, 572554, 572555, 573507, 573508, 573543, 573544, 573821, 573822, 575408, 575409, 575553, 575686, 575687, 577816, 577817, 578291, 578292, 580459, 580460, 581924, 581925, 582093, 582094, 583390, 583409, 589537, 591505, 591506, 593541, 593587, 593588, 593659, 593660, 594320, 594321, 594882, 596136, 596137, 598008, 598009, 598232, 598233, 598300, 598301, 600565, 600566, 602564, 603273, 603274, 603944, 603945, 604335, 604336, 606124, 606125, 606621, 606622, 607700, 607701, 609408, 609409, 613430, 613431, 613480, 613481, 617072, 617073, 622493, 622494, 624500, 624501, 624959, 624960, 626617, 628287, 628288, 628697, 630079, 633589, 633590, 633930, 633931, 635967, 636179, 636180, 637706, 637707, 637968, 641089, 641090, 641569, 641570, 643221, 643222, 643414, 643415, 646328, 646329, 648392, 650395, 650396, 653719, 653897, 653898, 654947, 654948, 655170, 655171, 656563, 656952, 656953, 657973, 657974, 658198, 658199, 658583, 658584, 659544, 659545, 659575, 659576, 659672, 659673, 662984, 667734, 667735, 670137, 670138, 673375, 673623, 673624, 674827, 674828, 675991, 675992, 676492, 676493, 677523, 677842, 678105, 679947, 679948, 680757, 680758, 681784, 681785, 685406, 685407, 687008, 687009, 688358, 688359, 688873, 688874, 689951, 689952, 690551, 691148, 691149, 692286, 692287, 693211, 693212, 694092, 694093, 695939, 695940, 698770, 700189, 700190, 700615, 700616, 700645, 700646, 700902, 703279, 703280, 705230, 705231, 705376, 705377, 706767, 706768, 706858, 706859, 707223, 707224, 707376, 707377, 713704, 713705, 715419, 715420, 715895, 717555, 721859, 721860, 723419, 725362, 727937, 728960, 728961, 730134, 731259, 733672, 733673, 735223, 735224, 736477, 736478, 737017, 737018, 738533, 738534, 739856, 739857, 740067, 740068, 741248, 741249, 741496, 741497, 743311, 743312, 745801, 747681, 747682, 747904, 747905, 752324, 752524, 752525, 752869, 759766, 759767, 761007, 761008, 762160, 762161, 762457, 762458, 764348, 764349, 766769, 766770, 767285, 767286, 768231, 768406, 768407, 768858, 768859, 771003, 771004, 773133, 774201, 774202, 775020, 775021, 776945, 776946, 777027, 777028, 777335, 777336, 777827, 777828, 778254, 778255, 778354, 778355, 778896, 778897, 781951, 781952,
        782272, 782273, 783539, 787206, 787207, 789461, 789462, 792312, 792313, 792699, 792700, 793931, 793932, 794036, 794037, 798078, 798079, 800722, 800723, 802218, 802219, 802353, 805849, 805850, 806788, 806789, 808031, 808032, 808385, 808386, 808387, 808388, 809201, 809202, 810497, 810498, 810843, 810895, 810896, 811748, 811749, 811823, 811824, 815779, 815780, 815941, 815942, 818977, 818978, 823772, 823773, 824065, 825523, 825524, 827692, 827693, 828499, 828500, 831493, 831494, 831731, 831732, 834252, 834253, 837257, 838819, 839126, 840504, 840505, 841781, 841782, 843764, 843822, 844524, 845377, 845378, 845473, 845474, 848154, 848155, 849136, 850390, 850391, 850502, 850503, 852663, 852664, 853522, 853523, 853925, 853926, 854833, 855270, 856728, 856729, 858158, 858159, 858522, 858523, 860057, 860058, 864407, 864408, 865265, 865266, 866769, 866770, 870109, 870110, 871503, 871504, 872383, 872929, 872930, 873891, 874363, 874364, 874692, 874693, 878214, 878215, 881335, 881336, 881458, 881459, 882352, 882353, 882607, 883436, 883437, 883522, 883523, 883847, 883848, 887755, 887756, 890157, 891444, 891747, 893694, 893695, 895424, 895425, 895726, 897107, 897108, 898693, 898694, 899246, 900818, 900819, 902132, 902133, 902996, 902997, 904534, 904535, 905083, 905084, 905245, 905246, 907850, 908644, 909319, 909320, 910791, 910792, 910793, 911120, 911121, 912924, 912925, 913438, 913439, 914160, 914161, 914454, 914455, 915604, 917253, 917254, 918939, 919257, 922316, 922317, 922408, 922409, 924778, 924779, 926292, 926293, 928656, 928657, 933064, 933065, 935921, 935922, 936083, 936084, 936536, 936537, 936775, 936776, 941277, 941278, 941836, 943814, 943815, 944531, 945214, 945215, 945586, 945587, 950094, 951485, 951486, 952301, 952302, 953033, 953034, 953115, 953116, 956618, 956619, 958280, 958281, 962221, 962222, 964577, 964578, 968596, 968597, 971966, 971967, 973575, 973576, 974508, 976426, 976427, 977236, 977818, 977819, 978061, 978062, 978359, 979921, 979922, 980895, 980896, 981658, 981659, 981925, 981926, 982268, 985749, 985750, 986355, 986356, 988471, 989017, 989018, 995269, 995270, 995570, 995571, 996248, 996249, 996615, 996616, 996870, 996871, 997172, 997173, 998480, 998481, 1000865, 1000866, 1001117, 1001118, 1002085, 1002086, 1002594, 1002595, 1003555, 1005370, 1005371, 1008073, 1008074, 1008196, 1008197, 1009293, 1009294, 1009754, 1009755, 1010269, 1010270, 1011112, 1011113, 1012336, 1015133, 1015134, 1017111, 1017112, 1018448, 1018449, 1019910, 1021776, 1021777, 1024018, 1024583, 1024584, 1026288, 1026298, 1026299, 1026460, 1026461, 1026880, 1026881, 1028211, 1028212, 1036000, 1036688, 1038550, 1038551, 1038879, 1038880, 1039806, 1039916, 1039917, 1039935, 1039936, 1040400,
        1040401, 1041107, 1041108, 1042981, 1044079, 1045606, 1045607, 1046256, 1047372, 1047373, 1047756, 1047757, 1048541, 1048542, 1050008, 1050009, 1053195, 1054241, 1054242, 1056203, 1056204, 1056649, 1056650, 1062183, 1062184, 1063419, 1065089, 1065090, 1065309, 1065310, 1065579, 1065580, 1067154, 1067155, 1072657, 1072658, 1077542, 1077543, 1078015, 1078016, 1080535, 1080536, 1080628, 1080629, 1080630, 1083517, 1083518, 1084621, 1084622, 1084673, 1084674, 1085109, 1085110, 1085335, 1085336, 1086106, 1086107, 1086517, 1087865, 1087866, 1088091, 1088092, 1088160, 1088673, 1088674, 1091707, 1091832, 1091833, 1091912, 1091913, 1092551, 1092552, 1093205, 1093206, 1096203, 1096204, 1096502, 1096503, 1097362, 1097363, 1097935, 1097936, 1098169, 1098170, 1099946, 1099947, 1100719, 1100720, 1100826, 1100827, 1100849, 1100850, 1101351, 1101352, 1104204, 1104205, 1104702, 1104703, 1105697, 1105698, 1105978, 1105979, 1112054, 1113360, 1113361, 1114646, 1114647, 1116668, 1116669, 1118339, 1118340, 1120075, 1120076, 1121379, 1121380, 1122708, 1122709, 1122888, 1122889, 1124339, 1124340, 1124969, 1124970, 1124984, 1124985, 1125426, 1129425, 1129426, 1132199, 1132692, 1132693, 1132833, 1132834, 1134129, 1134130, 1134507, 1134508, 1135824, 1137218, 1137219, 1137654, 1137655, 1138287, 1138288, 1139719, 1139720, 1142580, 1142581};
    std::sort(exclusionTable.begin(), exclusionTable.end());
    int entryRemovedNb = 0;
    int p0I = pSelListU.size() - 1;
    int p0J = exclusionTable.size() - 1;
    for(int i = p0I; i>0; i--){
        for(int j = p0J; j>-1; j--){
            // cout << exclusionTable[j]-pSelListU[i] << "\t" << pSelListU[i] <<":"<< exclusionTable[j]<< endl;
            if(pSelListU[i] > exclusionTable[j]){
                break;
            }else if(pSelListU[i] == exclusionTable[j]){
                // if(exclusionTable[j] == 1039) cout << "-\n";
                pSelListU.erase(pSelListU.begin() + i);
                entryRemovedNb++;
                p0J=j-1;
                break;
            }else{
                continue;
            }
        }        
    }
    // cout << endl << pSelListU.size() << endl;
    // for(int i=0; i<pSelListU.size(); i++){
    //     cout << pSelListU[i] << ", ";
    // }
    // cout << "end" << endl;
    
    // Downstream detector.
    // Selection happens using the track ID and the charge scoring plane tree
    class pSelListDCl {
        public:
        pSelListDCl(int event, int track){ eventID = event; trackID = track;};
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
            // cout << trackSelD << ", ";
            pSelListDCl tmp(eventSelD, trackSelD);
            pSelListD.push_back(tmp);
        }
    }
    Long64_t selDNb = pSelListD.size();
    if(verbosityLevel>0) cout << "Non-photons entering the downstream detector: " << selDNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListD.begin(), pSelListD.end(),
        // Lambda expression begins
        [](pSelListDCl a, pSelListDCl b) {
            return (a.getEvent() < b.getEvent() && a.getTrack() < b.getTrack());
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction
    // Debug
    // for(int i=0; i<selDNb; i++){
    //     cout << pSelListD[i].getEvent() << ":" <<  pSelListD[i].getTrack() << "\t";
    // }
    
    // Dose code with some minor modifications
    const int binsNb = 20./meshSize;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    TH2D* edepXYSelHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        edepXYSelHist[i] = new TH2D(TString::Format("edepXYSelHist%i",i), "non-photons primaries " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
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
                }else if(val == event){
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
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\t-chg\n";
                        edepXYSelHist[det]->Fill(x, y, edep);
                        pD0=j_events;
                        break;
                    }
                    continue;
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


    // Rescale to BX
    double factor = rescaleBX();
    if(factor == 0 && norm != 1){
        cout << "Manual normalization detected. This file corresponds to " << norm << " BX. Rescaling accordingly..." << endl;
        factor = norm * 1E-6;
    }
    
    //
    for(int i=0; i<2; i++){
        edepXYHist[i]->Scale(factor * 1E-6);
        edepXYSelHist[i]->Scale(factor * 1E-6);
        edepXYHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        edepXYSelHist[i]->GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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
        edepProj[i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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

    new TCanvas();
    double scale;
    scale = edepProj[4].GetXaxis()->GetBinWidth(1)/(edepProj[4].Integral());
    edepProj[4].Scale(scale);
    edepProj[4].DrawClone("hist");

    scale = edepProj[0].GetXaxis()->GetBinWidth(1)/(edepProj[0].Integral());
    edepProj[0].Scale(scale);
    edepProj[0].DrawClone("hist same");

    
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

    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        edepStrip[2*i].Scale(factor * 1E-6);
        edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        hitHist[i]->Scale(factor);
        hitHist[i]->GetYaxis()->SetTitle("hits");
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

    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        edepStrip[2*i].Scale(factor * 1E-6);
        edepStrip[2*i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
        hitHist[i].Scale(factor);
        hitHist[i].GetYaxis()->SetTitle("hits");
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
int createReport(TString path = "plots/xi5/", TString format = ".png", TString ofilename = "out.root"){
    int verbosityLevel_bak = verbosityLevel;
    verbosityLevel = -1;

    TFile* ofile = new TFile(ofilename, "recreate");
    TCanvas* painter = new TCanvas("painter", "Title");
    painter->SetGridx(0); painter->SetGridy(0);
    TString fname;
    TLegend* legend;
    
    std::vector<TH1D> temp;
    std::vector<TH2D> temp2;

    // PRIMARY BEAM
    // cout << "beamMonitor OK\n";
    // temp = beamMonitor();
    // for(int i=0; i < 3; i++){
    //     temp[i].Write();
    //     temp[i].Draw("hist");
    //     fname = temp[i].GetName();
    //     painter->Print(path+fname+format);
    // }
    // temp.clear();

    cout << "LUXE_beamMonitor OK\n";
    temp = beamMonitor();
    temp[0].SetLineColor(kRed); // X
    temp[1].SetLineColor(kBlue); // Y
    for(int i=1; i >= 0; i--){
        temp[i].Write();
        temp[i].Draw("hist same");
    }
    legend = new TLegend(0.75, 0.8, 0.90, 0.9);
    legend->SetBorderSize(0);   legend->SetFillColor(0);    legend->SetTextSize(0.04);
    legend->AddEntry(&temp[0], "X", "l");
    legend->AddEntry(&temp[1], "Y", "l");
    legend->Draw();
    fname = "LUXE_beamMonitor";
    painter->Print(path+fname+format);
    delete legend;
    temp.clear();


    // cout << "primaryProfile OK\n";
    // temp2 = primaryProfile();
    // temp2[0].Write();
    // temp2[0].Draw("LEGO1");
    // fname = temp2[0].GetName();
    // painter->Print(path+fname+format);
    // temp2.clear();
    

    // cout << "primaryEnergyProfile OK\n";
    // temp2 = primaryEnergyProfile();
    // temp2[0].Write();
    // temp2[0].Draw("LEGO2Z");
    // fname = temp2[0].GetName();
    // painter->Print(path+fname+format);
    // temp2.clear();
    

    // cout << "primarySpectrum OK\n";
    // temp = primarySpectrum(11, 16.0);
    // for(int i=0; i<2; i++){
    //     temp[i].Write();
    // }
    // TLegend* legend=new TLegend(0.75,0.8,0.90,0.9);
    // legend->SetBorderSize(0);   legend->SetFillColor(0);    legend->SetTextSize(0.04);
    // legend->AddEntry(&temp[1], "#gamma, others", "l");
    // legend->AddEntry(&temp[0], "e^{-}, e^{+}" , "l");
    // temp[0].SetStats(0); temp[1].SetStats(0);
    // temp[1].Draw(); temp[0].Draw("same");
    // legend->Draw();
    // fname = temp[0].GetName();
    // painter->Print(path+fname+format);
    // painter->SetLogx(); painter->SetLogy();
    // painter->Print(path+fname+"LogLog"+format);
    // painter->SetLogx(0); painter->SetLogy(0);
    // delete legend;
    // temp.clear();
    
    
    cout << "primarySpectrumLogLog OK\n";
    temp = primarySpectrumLogLog(1E-6, 20.0, 100);
    for(int i=0; i<2; i++){
        temp[i].Write();
    } 
    legend = new TLegend(0.25,0.8,0.45,0.9);
    legend->SetBorderSize(0);   legend->SetFillColor(0);    legend->SetTextSize(0.04);
    // legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
    legend->AddEntry(&temp[1], "#gamma, others", "l");
    legend->AddEntry(&temp[0], "others", "l");
    temp[0].SetStats(0); temp[1].SetStats(0);
    painter->SetLogx(); painter->SetLogy();
    temp[1].Draw("hist"); temp[0].Draw("hist same");
    legend->Draw();
    fname = temp[0].GetName();
    painter->Print(path+fname+format);
    painter->SetLogx(0); painter->SetLogy(0);
    delete legend;
    temp.clear();


    cout << "edepSpectrum OK\n";
    temp = edepSpectrum(500, 0, 500);
    for(int i=0; i<2; i++){
        temp[i].Write();
        temp[i].Draw();
        fname = temp[i].GetName();
        painter->Print(path+fname+format);
    }
    temp.clear();

    
    // These functions do not work properly. They must be fixed.
    // cout << "edepSpectrumChg OK\n";
    // temp = edepSpectrumChg(1000, 0, 1000);
    // TCanvas* edepSpectrumCanvas = new TCanvas("edepSpectrumCanvas", "Energy spectrum", 1200, 600);
    // edepSpectrumCanvas->DivideSquare(4);
    // for(int i=0; i<4; i++){
    //     edepSpectrumCanvas->cd(i+1);
    //     temp[i].Draw();
    // }
    // fname = "edepSpectrumChgMIP";
    // edepSpectrumCanvas->Print(path+fname+format);
    // delete edepSpectrumCanvas;
    // temp.clear();
    //

    cout << "edepStripChg OK\n";
    temp = edepStripChg();
    for(int i=0; i<2; i++){
        painter->SetName("edepStripChg"+whichDetShort(i));
        painter->SetTitle("Energy deposited per strip in the " + whichDetVerbose(i) + " detector");
        TLegend* legend=new TLegend(0.75,0.8,0.90,0.9);
        legend->SetBorderSize(0);   legend->SetFillColor(0);    legend->SetTextSize(0.04);
        legend->AddEntry(&temp[i], "#gamma, others", "l");
        legend->AddEntry(&temp[i+2], "others", "l");
        temp[i].SetStats(0); temp[i+2].SetStats(0);
        temp[i].Draw("hist");
        temp[i+2].Draw("hist same");
        legend->Draw();
        fname = painter->GetName();
        painter->Print(path+fname+format);
        delete legend;
    }
    temp.clear();
    // cout << "eSpectrumStrip OK\n";
    // temp = eSpectrumStrip(/*int stripNb = */100, /*float range = */60.);
    // for(int i=0; i<2; i++){
    //     temp[i].Write();
    //     temp[i].Draw();
    //     fname = temp[i].GetName();
    //     painter->Print(path+fname+format);
    // }
    // temp.clear();

    // cout << "eLongitudinal_array OK\n";
    // temp = eLongitudinal_array(/*nbBins = */100);
    // for(int i=0; i < 2; i++){
    //     temp[i].Write();
    //     temp[i].Draw("hist");
    //     fname = temp[i].GetName();
    //     painter->Print(path+fname+format);
    // }
    // temp.clear();

    here:
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
        temp2[i].SetStats(0);
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
    
    // cout << "edepProjectPrimEnPartitioned OK\n";
    // // Lambda function for formatting
    // auto formatEnergy = [](double x){
    //     TString unit;
    //     if(x>=0.1){
    //         unit = TString::Format("%.1f GeV", x);
    //     }else if(0.001<=x && x<0.1){
    //         unit = TString::Format("%.1f MeV", x*1E3);
    //     }else if(0.000001<=x && x<0.001){
    //         unit = TString::Format("%.1f keV", x*1E6);
    //     }else if(x<0.000001){
    //         unit = TString::Format("%.1f eV", x*1E9);
    //     }else{
    //         unit = TString::Format("%f", x);
    //     }
    //     return unit;
    // };
    // std::vector<TH1D> partitions = edepProjectPrimEnPartitioned();
    // int slices = partitions.size() / 2 - 1; //size is always even, so there is no possibility for int rounding errors
    // for(int det=0; det<2; det++){
    //     double maxCounts = 0;
    //     TLegend* edepProjectPrimEnPartitionedLegend = new TLegend(1.0-0.4, 0.7, 1.0-0.1, 0.9);
    //     edepProjectPrimEnPartitionedLegend->AddEntry(&partitions[2*slices+det], partitions[2*slices+det].GetTitle(), "l");//"0 < E < " + formatEnergy(partitions[2*slices+det].GetMaximum())
    //     for(int i= slices-1; i>=0; i--){
    //         edepProjectPrimEnPartitionedLegend->AddEntry(&partitions[2*i+det], partitions[2*i+det].GetTitle(), "l");
    //     }
    //     // Make sure that all histograms are viewed correctly
    //     for(int i= slices-1; i>=0; i--){
    //         if(partitions[2*i+det].GetMaximum() > maxCounts) maxCounts = partitions[2*i+det].GetMaximum();
    //     }
    //     partitions[2*slices+det].SetStats(0);
    //     partitions[2*slices+det].SetMinimum(1);
    //     partitions[2*slices+det].Draw("hist");
    //     for(int i= slices-1; i>=0; i--){
    //         partitions[2*i+det].GetYaxis()->SetRangeUser(0, 1.05*maxCounts); // Make sure that all histograms are viewed correctly
    //         partitions[2*i+det].SetStats(0);
    //         partitions[2*i+det].SetMinimum(1);
    //         partitions[2*i+det].Draw("hist same");
    //     }

    //     edepProjectPrimEnPartitionedLegend->Draw();
    //     fname = "edepProjectPrimEnPartitioned"+whichDetShort(det);
    //     painter->Print(path+fname+format);
    // }
    // partitions.clear();

    // // cout << "fitVSThrRng OK\n";
    // // std::vector<TGraph2D> fitVSThrRngPlots = fitVSThrRng();
    // // for(int i=0; i < 2; i++){
    // //     fitVSThrRngPlots[i].Write();
    // //     gStyle->SetPalette(1);
    // //     fitVSThrRngPlots[i].Draw("surf1");
    // //     fname = fitVSThrRngPlots[i].GetName();
    // //     painter->Print(path+fname+format);
    // // }
    // // fitVSThrRngPlots.clear();

    ofile->Close();
    verbosityLevel = verbosityLevel_bak;
    return 0;                                              
}


/*
* ******************************************************************************************************************************************
* NEW METHODS!!!
* ******************************************************************************************************************************************
*/
// int createReport(TString path = "plots/xi5/", TString format = ".png", TString ofilename = "out.root"){
//     int verbosityLevel_bak = verbosityLevel;
//     verbosityLevel = -1;

//     TFile* ofile = new TFile(ofilename, "recreate");
//     TCanvas* painter = new TCanvas("painter", "Title");
//     painter->SetGridx(); painter->SetGridy();
//     TString fname;

//     std::vector<TH1D> temp;
//     std::vector<TH2D> temp2;

//     // PRIMARY BEAM
//     // cout << "beamMonitor OK\n";
//     // temp = beamMonitor();
//     // for(int i=0; i < 3; i++){
//     //     temp[i].Write();
//     //     temp[i].Draw();
//     //     fname = temp[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // temp.clear();
    
//     // cout << "primaryProfile OK\n";
//     // temp2 = primaryProfile();
//     // temp2[0].Write();
//     // temp2[0].Draw("LEGO1");
//     // fname = temp2[0].GetName();
//     // painter->Print(path+fname+format);
//     // temp2.clear();
    
//     // cout << "primaryEnergyProfile OK\n";
//     // temp2 = primaryEnergyProfile();
//     // temp2[0].Write();
//     // temp2[0].Draw("LEGO2Z");
//     // fname = temp2[0].GetName();
//     // painter->Print(path+fname+format);
//     // temp2.clear();
    
//     // cout << "primarySpectrum OK\n";
//     // temp = primarySpectrum(11, 16.0);
//     // for(int i=0; i<2; i++){
//     //     temp[i].Write();
//     // }
//     // TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
//     // legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
//     // legend->AddEntry(&temp[1],TString::Format("#gamma,others- %i entries",(int)temp[1].GetEntries()),"l");
//     // legend->AddEntry(&temp[0],TString::Format("e^{-},e^{+}    - %i entries",(int)temp[0].GetEntries()),"l");
//     // temp[0].SetStats(0); temp[1].SetStats(0);
//     // temp[1].Draw(); temp[0].Draw("same");
//     // legend->Draw();
//     // fname = temp[0].GetName();
//     // painter->Print(path+fname+format);
//     // painter->SetLogx(); painter->SetLogy();
//     // painter->Print(path+fname+"LogLog"+format);
//     // painter->SetLogx(0); painter->SetLogy(0);
//     // delete legend;
//     // temp.clear();
    
//     cout << "primarySpectrumLogLog OK\n";
//     temp = primarySpectrumLogLog(11, 1E-9, 16.0, 100);
//     for(int i=0; i<2; i++){
//         temp[i].Write();
//     }
//     TLegend* legend4LogLog = new TLegend(0.1,0.7,0.4,0.9);
//     legend4LogLog->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
//     legend4LogLog->AddEntry(&temp[1],TString::Format("#gamma, others - %i entries",(int)temp[1].GetEntries()),"l");
//     legend4LogLog->AddEntry(&temp[0],TString::Format("others    - %i entries",(int)temp[0].GetEntries()),"l");
//     temp[0].SetStats(0); temp[1].SetStats(0);
//     painter->SetLogx(); painter->SetLogy();
//     temp[1].Draw(); temp[0].Draw("same");
//     legend4LogLog->Draw();
//     fname = temp[0].GetName();
//     painter->Print(path+fname+format);
//     painter->SetLogx(0); painter->SetLogy(0);
//     delete legend4LogLog;
//     temp.clear();


//     // cout << "edepSpectrum OK\n";
//     // temp = edepSpectrum(1000, 0, 1000);
//     // for(int i=0; i<2; i++){
//     //     temp[i].Write();
//     //     temp[i].Draw();
//     //     fname = temp[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // temp.clear();
    
    
//     // cout << "edepSpectrumChg OK\n";
//     // temp = edepSpectrumChg(1000, 0, 1000);
//     // TCanvas* edepSpectrumCanvas = new TCanvas("edepSpectrumCanvas", "Energy spectrum", 1200, 600);
//     // edepSpectrumCanvas->DivideSquare(4);
//     // for(int i=0; i<4; i++){
//     //     edepSpectrumCanvas->cd(i+1);
//     //     temp[i].Draw();
//     // }
//     // fname = "edepSpectrumChgMIP";
//     // edepSpectrumCanvas->Print(path+fname+format);
//     // delete edepSpectrumCanvas;
//     // temp.clear();
        

//     // cout << "edepStripChg OK\n";
//     // temp = edepStripChg(0);
//     // for(int i=0; i<2; i++){
//     //     painter->SetName("edepStripChg"+whichDetShort(i));
//     //     painter->SetTitle("Energy deposited per strip in the " + whichDetVerbose(i) + " detector");
//     //     TLegend* legend = new TLegend(1-0.4,0.7,1-0.1,0.9);
//     //     legend->SetHeader("Initial beam particles accounted:","L"); // option "C" allows to center the header
//     //     legend->AddEntry(&temp[i], TString::Format("#gamma, others - %i entries",(int)temp[i].GetEntries()),"l");
//     //     legend->AddEntry(&temp[i+2], TString::Format("others    - %i entries",(int)temp[i+2].GetEntries()),"l");
//     //     temp[i].SetStats(0); temp[i+2].SetStats(0);
//     //     temp[i].Draw("hist");
//     //     temp[i+2].Draw("hist same");
//     //     legend->Draw();
//     //     fname = painter->GetName();
//     //     painter->Print(path+fname+format);
//     //     delete legend;
//     // }
//     // temp.clear();


//     // cout << "eSpectrumStrip OK\n";
//     // temp = eSpectrumStrip(/*int stripNb = */100, /*float range = */60.);
//     // for(int i=0; i<2; i++){
//     //     temp[i].Write();
//     //     temp[i].Draw();
//     //     fname = temp[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // temp.clear();

//     // cout << "eLongitudinal_array OK\n";
//     // temp = eLongitudinal_array(/*nbBins = */100);
//     // for(int i=0; i < 2; i++){
//     //     temp[i].Write();
//     //     temp[i].Draw("hist");
//     //     fname = temp[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // temp.clear();

//     // cout << "edepProject OK\n";
//     // temp = edepProject(/*meshSize = */0.100);
//     // for(int i=0; i < 4; i++){
//     //     painter->cd(i);
//     //     temp[i].Write();
//     //     temp[i].Draw("hist");
//     //     fname = temp[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // temp.clear();

//     // cout << "doseXY OK\n";
//     // temp2 = doseXY();
//     // for(int i=0; i < 4; i++){
//     //     temp2[i].Write();
//     //     temp2[i].Draw("colz");
//     //     fname = temp2[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // temp2.clear();

//     // cout << "fitHits OK\n";
//     // std::vector<TH1I> fitHitsPlots = fitHits();
//     // for(int i=0; i < 2; i++){
//     //     fitHitsPlots[i].Write();
//     //     fitHitsPlots[i].Draw();
//     //     fname = fitHitsPlots[i].GetName();
//     //     painter->Print(path+fname+format);
//     // }
//     // fitHitsPlots.clear();
    
//     // cout << "edepProjectPrimEnPartitioned OK\n";
//     // // Lambda function for formatting
//     // auto formatEnergy = [](double x){
//     //     TString unit;
//     //     if(x>=0.1){
//     //         unit = TString::Format("%.1f GeV", x);
//     //     }else if(0.001<=x && x<0.1){
//     //         unit = TString::Format("%.1f MeV", x*1E3);
//     //     }else if(0.000001<=x && x<0.001){
//     //         unit = TString::Format("%.1f keV", x*1E6);
//     //     }else if(x<0.000001){
//     //         unit = TString::Format("%.1f eV", x*1E9);
//     //     }else{
//     //         unit = TString::Format("%f", x);
//     //     }
//     //     return unit;
//     // };
//     // std::vector<TH1D> partitions = edepProjectPrimEnPartitioned();
//     // int slices = partitions.size() / 2 - 1; //size is always even, so there is no possibility for int rounding errors
//     // for(int det=0; det<2; det++){
//     //     double maxCounts = 0;
//     //     TLegend* edepProjectPrimEnPartitionedLegend = new TLegend(1.0-0.4, 0.7, 1.0-0.1, 0.9);
//     //     edepProjectPrimEnPartitionedLegend->AddEntry(&partitions[2*slices+det], partitions[2*slices+det].GetTitle(), "l");//"0 < E < " + formatEnergy(partitions[2*slices+det].GetMaximum())
//     //     for(int i= slices-1; i>=0; i--){
//     //         edepProjectPrimEnPartitionedLegend->AddEntry(&partitions[2*i+det], partitions[2*i+det].GetTitle(), "l");
//     //     }
//     //     // Make sure that all histograms are viewed correctly
//     //     for(int i= slices-1; i>=0; i--){
//     //         if(partitions[2*i+det].GetMaximum() > maxCounts) maxCounts = partitions[2*i+det].GetMaximum();
//     //     }
//     //     partitions[2*slices+det].SetStats(0);
//     //     partitions[2*slices+det].SetMinimum(1);
//     //     partitions[2*slices+det].Draw("hist");
//     //     for(int i= slices-1; i>=0; i--){
//     //         partitions[2*i+det].GetYaxis()->SetRangeUser(0, 1.05*maxCounts); // Make sure that all histograms are viewed correctly
//     //         partitions[2*i+det].SetStats(0);
//     //         partitions[2*i+det].SetMinimum(1);
//     //         partitions[2*i+det].Draw("hist same");
//     //     }

//     //     edepProjectPrimEnPartitionedLegend->Draw();
//     //     fname = "edepProjectPrimEnPartitioned"+whichDetShort(det);
//     //     painter->Print(path+fname+format);
//     // }
//     // partitions.clear();

//     // // cout << "fitVSThrRng OK\n";
//     // // std::vector<TGraph2D> fitVSThrRngPlots = fitVSThrRng();
//     // // for(int i=0; i < 2; i++){
//     // //     fitVSThrRngPlots[i].Write();
//     // //     gStyle->SetPalette(1);
//     // //     fitVSThrRngPlots[i].Draw("surf1");
//     // //     fname = fitVSThrRngPlots[i].GetName();
//     // //     painter->Print(path+fname+format);
//     // // }
//     // // fitVSThrRngPlots.clear();

//     ofile->Close();
//     verbosityLevel = verbosityLevel_bak;
//     return 0;                                              
// }
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

// Energy deposition in the transverse XY plane. Either X or Y profile by selecting the returned value [0] det0X [1] det0Y [2] det1X [3] det1Y
std::vector<TH1D> ZZZZZ(bool filter=0, double meshSize = 0.100, double norm=1){
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
        if(pdgSelU != 22){
            pSelListU.push_back(eventSelU);
        }
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the upstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListU.begin(), pSelListU.end());

    
    
    // Downstream detector.
    // Selection happens using the track ID and the charge scoring plane tree
    class pSelListDCl {
        public:
        pSelListDCl(int event, int track){ eventID = event; trackID = track;};
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
            // cout << trackSelD << ", ";
            pSelListDCl tmp(eventSelD, trackSelD);
            pSelListD.push_back(tmp);
        }
    }
    Long64_t selDNb = pSelListD.size();
    if(verbosityLevel>0) cout << "Non-photons entering the downstream detector: " << selDNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListD.begin(), pSelListD.end(),
        // Lambda expression begins
        [](pSelListDCl a, pSelListDCl b) {
            return (a.getEvent() <= b.getEvent() && a.getTrack() < b.getTrack());
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction
    // Debug
    // for(int i=0; i<selDNb; i++){
    //     cout << pSelListD[i].getEvent() << ":" <<  pSelListD[i].getTrack() << "\t";
    // }
    // cout << endl;
    
    // Dose code with some minor modifications
    const int binsNb = 20./meshSize;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    TH2D* edepXYSelHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        edepXYSelHist[i] = new TH2D(TString::Format("edepXYSelHist%i",i), "non-photons primaries " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        // Style
        // edepXYHist[i]->SetStats(0);
        edepXYHist[i]->SetContour(500);
        edepXYSelHist[i]->SetContour(500);
        edepXYSelHist[i]->SetLineColor(kRed);
    }

    class debugTreeEntry{
        public:
        debugTreeEntry(int Revent, int Rtrack, int Rdet, int Rpdg, double Rx, double Ry, double Redep){
            event = Revent;
            track = Rtrack;
            det = Rdet;
            pdg = Rpdg;
            x = Rx;
            y = Ry;
            edep = Redep;
        };
        ~debugTreeEntry(){};

        int event;
        int track;
        int det;
        int pdg;
        double x;
        double y;
        double edep;
    };

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
    
    std::vector<debugTreeEntry> dTEntry;
    dTEntry.reserve(debugEntries);
    
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
        
        dTEntry.emplace_back(event, track, det, pdg, x, y, edep);
    }

    // cout << endl;
    // Debug
    // for(int i=0; i<dTEntry.size(); i++){
    //     cout << dTEntry[i].event << ":" <<  dTEntry[i].track << "\t";
    // }
    // cout << endl << "--------size: " << dTEntry.size() << " ----------debubugentries: "<< debugEntries << " -----------" << endl;


    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(dTEntry.begin(), dTEntry.end(),
        // Lambda expression begins
        [](debugTreeEntry a, debugTreeEntry b) {
            return (a.event < b.event);
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction

    // Understand which events are recorded
    std::vector<int> eventRecorded;
    int previousEvent = dTEntry[0].event;
    eventRecorded.emplace_back(0);
    for(int i=1; i<dTEntry.size(); i++){
        if(dTEntry[i].event != previousEvent){
            eventRecorded.push_back(i);
            previousEvent = dTEntry[i].event;
        }
    }
    eventRecorded.emplace_back(dTEntry.size());

    // Debug
    // for(int i=0; i< eventRecorded.size(); i++){
    //     cout << eventRecorded[i] << ", ";
    // }
    // cout << endl;

    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    for(int i=0; i< (eventRecorded.size()-1); i++){
        std::sort(dTEntry.begin() + eventRecorded[i], dTEntry.begin() + eventRecorded[i+1],
            // Lambda expression begins
            [](debugTreeEntry a, debugTreeEntry b) {
                return (a.track < b.track);
            } // end of lambda expression
        );    // The eventIDs are sorted by contruction
    }
    
    // Debug
    for(int i=0; i<dTEntry.size(); i++){
        if(dTEntry[i].det != 1) continue;
        cout << dTEntry[i].event << ":" <<  dTEntry[i].track << "\t";
    }
    cout << endl;
    
    cout << endl;
    cout << "----------------------------------------------------------------------------\n";
    // Debug
    for(int i=0; i<pSelListD.size(); i++){
        cout << pSelListD[i].getEvent() << ":" <<  pSelListD[i].getTrack() << "\t";
    }

    cout << "goodbye";
    // exit(-1);

    // Full the hist fol selected tracks
    for(int i=0; i<dTEntry.size(); i++){
        int event = dTEntry[i].event;
        int track = dTEntry[i].track;
        int det = dTEntry[i].det;
        int pdg = dTEntry[i].pdg;
        double x = dTEntry[i].x;
        double y = dTEntry[i].y;
        double edep = dTEntry[i].edep;
        // Filling energy depositions from selected 'filter' primaries

        cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\t";


        if(det==0){
            for(int j=pU0; j<selUNb; j++){
                int val=pSelListU[j];
                if(event < val){
                    break;
                }else if(val == event){
                    cout << "-chg";
                    edepXYSelHist[det]->Fill(x, y, edep);
                    // pU0=j;
                    break;
                }
            }
            cout << endl;
        }else{
            for(int j_events=pD0; j_events<pSelListD.size(); j_events++){
                auto pSelClass_ith = pSelListD[j_events];
                int evtOf_ith_class = pSelClass_ith.getEvent();
                int trkOf_ith_class = pSelClass_ith.getTrack();
                if(event < evtOf_ith_class){
                    break;
                }else if(event == evtOf_ith_class){
                    if(track == trkOf_ith_class){
                        cout << "-chg";
                        edepXYSelHist[det]->Fill(x, y, edep);
                        // pD0=j_events;
                        break;
                    }
                    continue;
                }
            }
            cout << endl;
        }
    }


    

    // Verbose messages
    if(verbosityLevel>1){
        cout << "Mesh size [um]: "  << meshSize*1000.           << endl;
        cout << "Mesh blocks: "     << binsNb<<'x'<<binsNb      << endl;
        cout    << "\n\n";
    }


    // Rescale to BX
    double factor = rescaleBX();
    for(int i=0; i<2; i++){
        edepXYHist[i]->Scale(factor * 1E-6);
        edepXYSelHist[i]->Scale(factor * 1E-6);
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
        edepProj[i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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





// Export to file the particles crossing the chg. sensitive plane
void dumpParticlesCrossingDownsteramGBP(){
    TFile* newfile = new TFile("kyleBeam_6_50_0_w100_downGBP.root", "RECREATE");
	// Create new TTree to fill
    Int_t pdg;
	Double_t energy, weight, xx, yy, zz, pxx, pyy, pzz;
	TTree* SimTrack = new TTree("Tracks", "MCTracks crossing downstream GBP");
	SimTrack->Branch("pdg", &pdg, "pdg/I");
	SimTrack->Branch("x", &xx, "x/D");
	SimTrack->Branch("y", &yy, "y/D");
	SimTrack->Branch("z", &zz, "z/D");
	SimTrack->Branch("px", &pxx, "px/D");
	SimTrack->Branch("py", &pyy, "py/D");
	SimTrack->Branch("pz", &pzz, "pz/D");
	SimTrack->Branch("energy", &energy, "energy/D");
	SimTrack->Branch("weight", &weight, "weight/D");

    int w0 = -getConfig();
    int eventSelD, trackSelD, pdgSelD;
    double x0, y0, z0, px, py, pz, ekin;
    planeTree->SetBranchAddress("pdg", &pdgSelD);
    planeTree->SetBranchAddress("x0", &x0);
    planeTree->SetBranchAddress("y0", &y0);
    planeTree->SetBranchAddress("z0", &z0);
    planeTree->SetBranchAddress("px", &px);
    planeTree->SetBranchAddress("py", &py);
    planeTree->SetBranchAddress("pz", &pz);
    planeTree->SetBranchAddress("ekin", &ekin);
    Long64_t planeEntries = planeTree->GetEntries();
    if(verbosityLevel>0) cout << "planeTree has " << planeEntries << " entries. Looping over entries...\n";
    resetProgress();
    for(Long64_t i=0; i< planeEntries; i++ ){
        float status = (float)(i+1) / planeEntries;
        printProgress(status);
        planeTree->GetEntry(i);

        pdg    = pdgSelD;
		xx     = x0;
		yy     = y0;
		zz     = z0; // mm
		pxx    = px;
		pyy    = py;
		pzz    = pz;
		energy = ekin;
		weight = w0;
		SimTrack->Fill();
    }
    SimTrack->Write();
	delete SimTrack;
	newfile->Close();
}


// This is the correct function for accounting charged particles downstream
// Energy deposition in the transverse XY plane. Either X or Y profile by selecting the returned value [0] det0X [1] det0Y [2] det1X [3] det1Y
std::vector<TH1D> XXXX_LIFEBOAT(bool filter=0, double meshSize = 0.100, double norm=1){
    // SELECTION FOR THE UPSTREAM DETECTOR - USING THE EVENT TREE BECAUSE PARTICLES ARE SPAWNED THERE
    // SELECTION FOR THE UPSTREAM DETECTOR - USING THE EVENT TREE BECAUSE PARTICLES ARE SPAWNED THERE
    // SELECTION FOR THE UPSTREAM DETECTOR - USING THE EVENT TREE BECAUSE PARTICLES ARE SPAWNED THERE
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
        if(pdgSelU != 22){
            pSelListU.push_back(eventSelU);
        }
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the downstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListU.begin(), pSelListU.end());


// SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE// SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE// SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE// SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE

    // SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE
    // SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE
    // SELECTION FOR THE DOWNSTREAM DETECTOR - USING THE PLANE TREE
    // Selection happens using the track ID and the charge scoring plane tree
    class pSelListDCl {
        public:
        pSelListDCl(int revent, int rtrack){ event = revent; track = rtrack;};
        ~pSelListDCl(){};

        int getEvent(){ return event;}
        int getTrack(){ return track;}

        int event;
        int track;
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
            // cout << trackSelD << ", ";
            pSelListDCl tmp(eventSelD, trackSelD);
            pSelListD.push_back(tmp);
        }
    }
    Long64_t selDNb = pSelListD.size();
    if(verbosityLevel>0) cout << "Non-photons entering the downstream detector: " << selDNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListD.begin(), pSelListD.end(),
        // Lambda expression begins
        [](pSelListDCl a, pSelListDCl b) {
            return (a.getEvent() < b.getEvent());
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction
    std::vector<int> pSelListDeventRecorded;
    int pSelListpreviousEvent = pSelListD[0].getEvent();
    pSelListDeventRecorded.emplace_back(0);
    for(int i=1; i<pSelListD.size(); i++){
        if(pSelListD[i].getEvent() != pSelListpreviousEvent){
            pSelListDeventRecorded.push_back(i);
            pSelListpreviousEvent = pSelListD[i].getEvent();
        }
    }
    pSelListDeventRecorded.emplace_back(pSelListD.size()); // dTEntry.size() is the last element of the array
    // Sort the list of SELECTED-aka filtered part IDs in order to speed up the loop by using p0=i
    for(int i=0; i< (pSelListDeventRecorded.size()-1); i++){
        std::sort(pSelListD.begin() + pSelListDeventRecorded[i], pSelListD.begin() + pSelListDeventRecorded[i+1],
            // Lambda expression begins
            [](pSelListDCl a, pSelListDCl b) {
                return (a.getTrack() < b.getTrack());
            } // end of lambda expression
        );    // The eventIDs are sorted by contruction
    }
    // Debug
    // for(int i=0; i<selDNb; i++){
    //     cout << pSelListD[i].getEvent() << ":" <<  pSelListD[i].getTrack() << "\t";
    // }
    // cout << endl;

    
    // Dose code with some minor modifications
    const int binsNb = 20./meshSize;
    double totEnergy[2] = {0};

    TH2D* edepXYHist[2];
    TH2D* edepXYSelHist[2];
    for(int i=0; i<2; i++){
        edepXYHist[i] = new TH2D(TString::Format("edepXYHist%i",i), "Energy dep. XY map in the " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        edepXYSelHist[i] = new TH2D(TString::Format("edepXYSelHist%i",i), "non-photons primaries " + whichDetVerbose(i) + " detector;X [mm]; Y [mm]", binsNb, -10, 10, binsNb, -10, 10);
        // Style
        // edepXYHist[i]->SetStats(0);
        // edepXYHist[i]->SetContour(500);
        // edepXYSelHist[i]->SetContour(500);
        edepXYSelHist[i]->SetLineColor(kRed);
    }

    class debugTreeEntry{
        public:
        debugTreeEntry(int Revent, int Rtrack, int Rptrack, int Rdet, int Rpdg, double Rx, double Ry, double Redep){
            event = Revent;
            track = Rtrack;
            ptrack = Rptrack;
            det = Rdet;
            pdg = Rpdg;
            x = Rx;
            y = Ry;
            edep = Redep;
        };
        ~debugTreeEntry(){};

        int event;
        int track;
        int ptrack;
        int det;
        int pdg;
        double x;
        double y;
        double edep;
    };

    if(verbosityLevel>0) cout << "Filling the histogram..." << endl;
    int event, track, ptrack, det, pdg; double x, y; double edep;
    debugTree->SetBranchAddress("event", &event);           // Event ID
    debugTree->SetBranchAddress("track", &track);           // Track ID
    debugTree->SetBranchAddress("ptrack", &ptrack);         // Track ID
    debugTree->SetBranchAddress("det", &det);               // 0,1
    debugTree->SetBranchAddress("x", &x);                   // mm
    debugTree->SetBranchAddress("y", &y);                   // mm
    debugTree->SetBranchAddress("edep", &edep);             // keV
    debugTree->SetBranchAddress("pdg", &pdg);               // -11, 11, 22, -
    Long64_t debugEntries = debugTree->GetEntries();
    
    std::vector<debugTreeEntry> dTEntry;
    dTEntry.reserve(debugEntries);
    
    // First loop, create the energy deposition map
    // AND save each deposition in the dTEntry VECTOR
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
        
        dTEntry.emplace_back(event, track, ptrack, det, pdg, x, y, edep);        
    }

    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(dTEntry.begin(), dTEntry.end(),
        // Lambda expression begins
        [](debugTreeEntry a, debugTreeEntry b) {
            return (a.event < b.event);
        } // end of lambda expression
    );    // The eventIDs are sorted by contruction

    // Create a list with the entry number of the degubTree where an edep from a different event occurred.
    // per esempio le depositioni associatea ll'evento 0 vanno dall'entrata 0 alla 3, quelle dell'evento 1 dalla 4 alla 9 eccetera
    // questa funzione salva 0, 3, 9 eccetera... fino alla fine debugTreeEntries !
    // questa cosa serve per fare poi il sort delle deposizioni evento per evento basandosi sul numero di traccia
    std::vector<int> eventRecorded;
    int previousEvent = dTEntry[0].event;
    eventRecorded.emplace_back(0);
    for(int i=1; i<dTEntry.size(); i++){
        if(dTEntry[i].event != previousEvent){
            eventRecorded.push_back(i);
            previousEvent = dTEntry[i].event;
        }
    }
    eventRecorded.emplace_back(dTEntry.size()); // dTEntry.size() is the last element of the array
    // Event by event,
    // Sort the list of SELECTED-aka filtered part IDs in order to speed up the loop by using p0=i
    for(int i=0; i< (eventRecorded.size()-1); i++){
        std::sort(dTEntry.begin() + eventRecorded[i], dTEntry.begin() + eventRecorded[i+1],
            // Lambda expression begins
            [](debugTreeEntry a, debugTreeEntry b) {
                return (a.track < b.track);
            } // end of lambda expression
        );    // The eventIDs are sorted by contruction
    }
    




    //Update the tables for positions of events (because the sorting changed positions)
    eventRecorded.clear();
    previousEvent = dTEntry[0].event;
    eventRecorded.emplace_back(0);
    for(int i=1; i<dTEntry.size(); i++){
        if(dTEntry[i].event != previousEvent){
            eventRecorded.push_back(i);
            previousEvent = dTEntry[i].event;
        }
    }
    eventRecorded.emplace_back(dTEntry.size());
    //
    pSelListDeventRecorded.clear();
    pSelListpreviousEvent = pSelListD[0].getEvent();
    pSelListDeventRecorded.emplace_back(0);
    for(int i=1; i<pSelListD.size(); i++){
        if(pSelListD[i].getEvent() != pSelListpreviousEvent){
            pSelListDeventRecorded.push_back(i);
            pSelListpreviousEvent = pSelListD[i].getEvent();
        }
    }
    pSelListDeventRecorded.emplace_back(pSelListD.size()); // dTEntry.size() is the last element of the array


    // first, extract from the edep list the edeps from secondary particles
    std::vector<debugTreeEntry> secondary_dTEntry;
    for(int i=0; i<dTEntry.size(); i++){
        if(dTEntry[i].ptrack > 0){
            secondary_dTEntry.push_back(dTEntry[i]);
        }
    }




    //Understand where events finish and order by ptrackid
    std::vector<int> eventRecordedC;
    previousEvent = secondary_dTEntry[0].event;
    eventRecordedC.emplace_back(0);
    for(int i=1; i<secondary_dTEntry.size(); i++){
        if(secondary_dTEntry[i].event != previousEvent){
            eventRecordedC.push_back(i);
            previousEvent = secondary_dTEntry[i].event;
        }
    }
    eventRecordedC.emplace_back(secondary_dTEntry.size()); // secondary_dTEntry.size() is the last element of the array
    // Event by event,
    // Sort the list of SELECTED-aka filtered part IDs in order to speed up the loop by using p0=i
    for(int i=0; i< (eventRecordedC.size()-1); i++){
        std::sort(secondary_dTEntry.begin() + eventRecordedC[i], secondary_dTEntry.begin() + eventRecordedC[i+1],
            // Lambda expression begins
            [](debugTreeEntry a, debugTreeEntry b) {
                return (a.ptrack < b.ptrack);
            } // end of lambda expression
        );    // The eventIDs are sorted by contruction
    }
    // Debug
    for(int i=0; i< secondary_dTEntry.size()-1; i++){
        cout << secondary_dTEntry[i].event << ":" <<  secondary_dTEntry[i].track << ":" << secondary_dTEntry[i].ptrack << "\t";
    }
    cout << endl;
    cout << "Secondary tracks: " << secondary_dTEntry.size() << endl;

    class toEraseClass{
        public:
            toEraseClass(debugTreeEntry rJ, debugTreeEntry rK, int rMoth):
            dad(rJ),
            son(rK),
            grandma(rMoth) {};
            ~toEraseClass(){};

            debugTreeEntry dad;
            debugTreeEntry son;
            int grandma;
    };

    std::vector<toEraseClass> toErase;
    // Event by event, trace back the particle to its mother. Then compare mother track ID with the one recorded in the plane
    for(int i=0; i<eventRecordedC.size()-1; i++){
        for(int j=eventRecordedC[i+1]-1; j>eventRecordedC[i]; j--){
            int previousJ = j;
            int k=-1;
            int grandma_min = secondary_dTEntry[j].ptrack;
            while(j+k > 0){
                if(secondary_dTEntry[j].ptrack == secondary_dTEntry[j+k].track){
                    // cout << j << "\t" << k << "\t" << secondary_dTEntry[j].track<<":"<<secondary_dTEntry[j].ptrack << " " << secondary_dTEntry[j+k].track<<":"<<secondary_dTEntry[j+k].ptrack << endl;
                    toErase.emplace_back(secondary_dTEntry[j], secondary_dTEntry[j+k], std::min(secondary_dTEntry[j+k].ptrack, grandma_min));
                    j=j+k;
                }else{
                    k--;
                }
            }
            // if(toErase.size()>0){
            //     const toEraseClass& grandMother = toErase.back();
            //     cout << "Event: " << i << " " << grandMother.mother << endl;
            //     toErase.clear();
            // }
            j=previousJ;
        }
    }

    for(int i=0; i<toErase.size(); i++){
        cout << "Event- " << toErase[i].dad.event<<":"<<toErase[i].dad.track<<":"<<toErase[i].dad.ptrack << "  " << toErase[i].son.event<<":"<<toErase[i].son.track<<":"<<toErase[i].son.ptrack << "  grandma: " << toErase[i].grandma << endl;
    }
    exit(0);
    // // Trace back the secondary to its grandmother
    // for(int i=0; i<secondary_dTEntry.size(); i++){
    //     for(int j=0; j<pSelListD.size(); j++){
    //         pSelListD[j].track
    //     }
    // }



    








   

    // Debug
    // for(int i=0; i<dTEntry.size(); i++){
    //     if(dTEntry[i].det != 1) continue;
    //     cout << dTEntry[i].event << ":" <<  dTEntry[i].track << ":" << dTEntry[i].ptrack << "\t";
    //     // cout << dTEntry[i].event << ":" <<  dTEntry[i].track << "\t";
    // }
    // cout << endl;
    
    // cout << endl;
    // cout << "----------------------------------------------------------------------------\n";
    // // Debug
    // for(int i=0; i<pSelListD.size(); i++){
    //     cout << pSelListD[i].getEvent() << ":" <<  pSelListD[i].getTrack() << "\t";
    // }



    cout << "goodbye";
    exit(-1);

    // Full the hist fol selected tracks - DET 1
    for(int i=0; i<dTEntry.size(); i++){
        int event = dTEntry[i].event;
        int track = dTEntry[i].track;
        int det = dTEntry[i].det;
        int pdg = dTEntry[i].pdg;
        double x = dTEntry[i].x;
        double y = dTEntry[i].y;
        double edep = dTEntry[i].edep;
        // Filling energy depositions from selected 'filter' primaries

        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\t";
        if(det==1){
            for(int j=pU0; j<selUNb; j++){
                int val=pSelListU[j];
                if(event < val){
                    // break;
                }else if(val == event){
                    // cout << "-chg";
                    edepXYSelHist[det]->Fill(x, y, edep);
                    // pU0=j;
                    // break;
                }
            }
            // cout << endl;
        }
    }

    // Rescale to BX
    // double factor = rescaleBX();
    // for(int i=0; i<2; i++){
    //     edepXYHist[i]->Scale(factor * 1E-6);
    //     edepXYSelHist[i]->Scale(factor * 1E-6);
    // }


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
        edepProj[i].GetYaxis()->SetTitle("Deposited energy / BX [GeV]");
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

// Same but for the spectrum
std::vector<TH1D> YYYYYY_edepSpectrumChg_DOWNSTREAM(int bins = 500, double a=0, double b=500){
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
        if(pdgSelU != 22){
            pSelListU.push_back(eventSelU);
        }
    }
    Long64_t selUNb = pSelListU.size();
    if(verbosityLevel>0) cout << "Non-photons entering the upstream detector: " << selUNb <<  endl;
    // Sort the list of chg. part IDs in order to speed up the loop by using p0=i
    std::sort(pSelListU.begin(), pSelListU.end());


    // Dose code with some minor modifications
    double totEnergy[2] = {0};

    TH1D* eSpectrHist[2];
    TH1D* eSpectrSelHist[2];
    for(int i=0; i<2; i++){
        eSpectrHist[i] = new TH1D(TString::Format("edepSpectrDet%i",i), "Energy spectrum " + whichDetVerbose(i) + " detector;keV; counts", bins, a, b);
        eSpectrSelHist[i] = new TH1D(TString::Format("edepSpectrSelDet%i",i), "Energy spectrum " + whichDetVerbose(i) + " detector pdg!=22;keV; counts", bins, a, b);
        // Style
        eSpectrSelHist[i]->SetLineColor(kRed);
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

    int prevEvent=1;
    double edepInEvent[2] = {};
    double edepInEventSel[2] = {};

    int pU0=0; int pD0=0;
    resetProgress();
    for(Long64_t i=0; i< debugEntries; i++ ){
        debugTree->GetEntry(i);
        float status = (float)(i+1)/(float)debugEntries;
        printProgress(status);

        if(det >= 2) continue; // Ignore the charge scorer plane
        totEnergy[det] += edep;
        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " << edep << "\t";

        // The event is the same of the previous one. Keep summing energy dep. contribution in the same arrays (edepInEvent or edepInEventSel)
        if(event == prevEvent){
            // cout << "=";
            // Energy contributions from all tracks without any tagging
            edepInEvent[det] += edep;

            // Energy contribution from selected tracks
            if(det==1){
                // cout << " det0\n";
                // For det0 selection is made using eventID since only 1 track is injected in sapphire because of the way the beam is imported
                for(int j=pU0; j<selUNb; j++){
                    int val=pSelListU[j];
                    if(val > event){
                        // break;
                    }else if(val == event){
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\tchg\n";
                        edepInEventSel[det] += edep;
                        // pU0=j;
                        // break;
                    }
                }
            }
        }else{
            // cout << "- ";
            // The event changed. Store the value for the edep and continue
            for(int ii=0; ii<2; ii++){
                // cout << edepInEvent[ii] << ", " << edepInEventSel[ii] << "; ";
                if(edepInEvent[ii] != 0){
                    eSpectrHist[ii]->Fill(edepInEvent[ii]);
                    // cout << "filling htot with: " << edepInEvent[ii] << endl;
                    edepInEvent[ii] = 0;
                }
                if(edepInEventSel[ii] != 0){
                    eSpectrSelHist[ii]->Fill(edepInEventSel[ii]);
                    // cout << "filling hSel with: " << edepInEventSel[ii] << endl;
                    edepInEventSel[ii] = 0;
                }
            }

            // Add this entry deposition
            edepInEvent[det] += edep;
            // Energy contribution from selected tracks
            if(det==1){
                // cout << " det0\n";
                // For det0 selection is made using eventID since only 1 track is injected in sapphire because of the way the beam is imported
                for(int j=pU0; j<selUNb; j++){
                    int val=pSelListU[j];
                    if(val > event){
                        // break;
                    }else if(val == event){
                        // cout << i << ": " << event << ", " << track << ", " << whichDetShort(det) << ", " << pdg << ", (" << x << "," << y << "), " <<  edep << "\tchg\n";
                        edepInEventSel[det] += edep;
                        // pU0=j;
                        // break;
                    }
                }
            }
            prevEvent = event;
        }
    }

    if(verbosityLevel>0){
        cout << "tot edep in det0: " << totEnergy[1] << endl; 
        TCanvas* edepSpectrumCanvas = new TCanvas("edepSpectrumCanvas", "Energy spectrum", 1200, 600);
        edepSpectrumCanvas->Divide(2,2);
        edepSpectrumCanvas->cd(1);
        eSpectrHist[0]->Draw();
        edepSpectrumCanvas->cd(2);
        eSpectrHist[1]->Draw();
        edepSpectrumCanvas->cd(3);
        eSpectrSelHist[0]->Draw();
        edepSpectrumCanvas->cd(4);
        eSpectrSelHist[1]->Draw();
    }

    std::vector<TH1D> retVar = {TH1D(*eSpectrHist[0]), TH1D(*eSpectrHist[1]), TH1D(*eSpectrSelHist[0]), TH1D(*eSpectrSelHist[1]) };
    return retVar;
}

// Superimpose the beam profile of the primaries with the reconstructed one
void SIGOVERLAY(){
    TFile* signalFile = new TFile("import/geant_src/e0lp_5_0_0_particles_g4_kyleTracks_signal.root", "read");
    TTree* signalTree = (TTree*)signalFile->Get("Tracks");
    TH1D* signalXHist = new TH1D("signalXHist", "signalXHist", 200, -10, 10);
    TH1D* signalYHist = new TH1D("signalYHist", "signalYHist", 200, -10, 10);
    signalTree->Project("signalXHist", "x");
    signalTree->Project("signalYHist", "y");
    signalXHist->SetLineStyle(2);
    signalXHist->SetLineColorAlpha(kRed, 0.5);
    signalYHist->SetLineStyle(2);
    signalYHist->SetLineColorAlpha(kRed, 0.5);


    cout << "edepProject OK\n";
    std::vector<TH1D> temp = edepProject(); // meshSize=0.100
    temp[0].SetLineColorAlpha(kBlack, 0.5);
    temp[3].SetLineColorAlpha(kBlack, 0.5);

    signalXHist->Scale(signalXHist->GetXaxis()->GetBinWidth(1)/(signalXHist->Integral()));
    temp[0].Scale(temp[0].GetXaxis()->GetBinWidth(1)/(temp[0].Integral()));

    signalYHist->Scale(signalYHist->GetXaxis()->GetBinWidth(1)/(signalYHist->Integral()));
    temp[3].Scale(temp[3].GetXaxis()->GetBinWidth(1)/(temp[3].Integral()));

    double maxX = std::max(signalXHist->GetMaximum(), temp[0].GetMaximum());
    signalXHist->GetYaxis()->SetRangeUser(0, 1.05*maxX);
    temp[0].GetYaxis()->SetRangeUser(0, 1.05*maxX);

    double maxY = std::max(signalYHist->GetMaximum(), temp[3].GetMaximum());
    signalYHist->GetYaxis()->SetRangeUser(0, 1.05*maxY);
    temp[3].GetYaxis()->SetRangeUser(0, 1.05*maxY);

    signalXHist->GetYaxis()->SetLabelSize(0);
    temp[0].GetYaxis()->SetLabelSize(0);
    signalYHist->GetYaxis()->SetLabelSize(0);
    temp[3].GetYaxis()->SetLabelSize(0);

    TCanvas* SIGOVERLAYCanvas = new TCanvas("SIGOVERLAYCanvas", "Overlay of measured profile vs compton one", 1200, 400);
    SIGOVERLAYCanvas->Divide(2);
    SIGOVERLAYCanvas->cd(1);
    TLegend* uLegend = new TLegend(0.75, 0.7, 0.90, 0.9);
    uLegend->SetBorderSize(0);   uLegend->SetFillColor(0);    uLegend->SetTextSize(0.04);
    uLegend->SetHeader("upstream","L"); // option "C" allows to center the header
    uLegend->AddEntry(signalXHist, "#gamma compton", "l");
    uLegend->AddEntry(&temp[0], "profiler", "l");
    signalXHist->Draw("hist");
    temp[0].DrawClone("hist same");
    uLegend->Draw();

    SIGOVERLAYCanvas->cd(2);
    TLegend* dLegend = new TLegend(0.75, 0.7, 0.90, 0.9);
    dLegend->SetBorderSize(0);   dLegend->SetFillColor(0);    dLegend->SetTextSize(0.04);
    dLegend->SetHeader("downstream","L"); // option "C" allows to center the header
    dLegend->AddEntry(signalXHist, "#gamma compton", "l");
    dLegend->AddEntry(&temp[3], "profiler", "l");
    signalYHist->Draw("hist");
    temp[3].DrawClone("hist same");
    dLegend->Draw();
}