void readTVMAdata_fit(TString inFile = "tmva_graph_data.root")
{
    TFile* tf = new TFile(inFile);
    TTree *traintr = (TTree*)tf->Get("trainTree");

    Float_t lnuEnergy = 0.;
    Float_t muon_mom = 0;
    Float_t totPh = 0.;
    Float_t totCubes = 0;

    traintr->SetBranchAddress("totalCubes", &totCubes);
    traintr->SetBranchAddress("totalPhotons", &totPh);
    traintr->SetBranchAddress("nuEnergy", &lnuEnergy);
    traintr->SetBranchAddress("muon_mom", &muon_mom);

    TTree *fittr = (TTree*)tf->Get("fittedMomTree");
    
    Float_t fit_lnuEnergy = 0.;
    Float_t fit_muon_mom = 0;
    Float_t fit_totPh = 0.;
    Float_t fit_totCubes = 0;

    fittr->SetBranchAddress("totalCubes", &fit_totCubes);
    fittr->SetBranchAddress("totalPhotons", &fit_totPh);
    fittr->SetBranchAddress("nuEnergy", &fit_lnuEnergy);
    fittr->SetBranchAddress("muon_mom", &fit_muon_mom);



    int entries = traintr->GetEntries();
    double sum(0.);
    int count(0);
    for(int i=0; i < entries; ++i) 
    { 
        //cout << "               Event " << i << "           "<< endl;
        traintr->GetEntry(i); 
        fittr->GetEntry(i);

        if(muon_mom!=0)
        {
            double percentage = std::fabs(100. *(muon_mom - fit_muon_mom)/ muon_mom);
            sum += percentage;
            count++;
            cout << " MC muon mom " << muon_mom << " , Fitted Mom " << fit_muon_mom  << "[ " << percentage << " per]"<< endl;
        }
        
    }
    cout << " ======================================="<< endl;
    cout << " Avg percentage " << sum/count << endl;
}
