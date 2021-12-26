// Progress monitor
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
float globalProgress = 0; bool globalProgressFlag = true;
void resetProgress(){ globalProgress = 0; globalProgressFlag=true; }
void printProgress(float percentage) {
    if( (percentage-globalProgress) >= 0.01) {
        int val = (int) (percentage * 100 + 1);
        int lpad = (int) (percentage * PBWIDTH + 1);
        int rpad = PBWIDTH - lpad;
        printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
        fflush(stdout);
        globalProgress = percentage;
    }else if(percentage >= 1. && globalProgressFlag){
        globalProgressFlag = false;
        cout << endl;
    }
}

// Export an histogram in a text file. Default separator is tabulation
void exportDAT(TH1D* hist, TString title = "filename.txt", TString separator = "\t"){
    ofstream myfile;
    myfile.open(title);
    myfile << "bincenter" << separator << "value" << endl;
    // How many bins
    Long64_t bins = hist->GetNbinsX();
    
    for(Long64_t index=0; index < bins; index++){
        double x = hist->GetBinCenter(index);
        double val = hist->GetBinContent(index);
        myfile << x << separator << val;
        if(index!=(bins-1)) myfile << endl;
    }
    myfile.close();
}

// Export an histogram in a text file using CSV format
void exportCSV(TH1D* hist, TString title = "filename.txt"){
    exportDAT(hist, title, ",");
}

// Guess the configuration from the filename
int getConfig(TString name){
    if(name.Contains("_9_")){
        return 9;
    }else if(name.Contains("_7_")){
        return 7;
    }else if(name.Contains("_6_") || name.Contains("_sbr_")){
        return 6;
    }else{
        return -1;
    }
}