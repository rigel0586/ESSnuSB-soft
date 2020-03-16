void tmva_2_total_ph_app(TString inFile = "tmva_data.root")
{
    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 

    Float_t totalPhotons, totalCubes;
    reader->AddVariable( "totalPhotons", &totalPhotons );
    reader->AddVariable( "totalCubes", &totalCubes );

    TString dir    = "dataset/weights/";
    TString prefix = "TMVARegression";
    TString methodName = "MLP";
    TString weightfile = dir + prefix + "_" + methodName + ".weights.xml";
    reader->BookMVA( methodName, weightfile ); 

    TFile* tf = new TFile(inFile);
    TTree *tr = (TTree*)tf->Get("TotalPhotonsTree");

    tr->SetBranchAddress( "totalPhotons", &totalPhotons );
    tr->SetBranchAddress( "totalCubes", &totalCubes );

    

    TCanvas* canvas = new TCanvas();
    TH1F* hist = new TH1F("nu E", "Neutrino Energy", 200, 0, 2);

    Int_t limit = tr->GetEntries();
    for(Int_t j = 0 ; j < limit; ++j)
    {
        tr->GetEntry(j); 
        std::vector<Float_t> result = reader->EvaluateRegression( methodName );
        for(Int_t i = 0; i < result.size(); ++i)
        {
            if(result[i]!=0)
                hist->Fill( result[i] );
        }
    }

    hist->Draw("colz");
    canvas->SaveAs("nu_E");

	delete reader;

	return 0;
}
