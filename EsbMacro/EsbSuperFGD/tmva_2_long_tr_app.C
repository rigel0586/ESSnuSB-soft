void tmva_2_long_tr_app(TString inFile = "tmva_data.root", TString methodName = "MLP")
{
    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 

    Float_t tr1, tr2, tr3, ph_tr1, ph_tr2, ph_tr3;
    reader->AddVariable( "tr1", &tr1 );
    reader->AddVariable( "tr2", &tr2 );
    reader->AddVariable( "tr3", &tr3 );
    reader->AddVariable( "ph_tr1", &ph_tr1 );
    reader->AddVariable( "ph_tr2", &ph_tr2 );
    reader->AddVariable( "ph_tr3", &ph_tr3 );

    Float_t totalPhotons, totalCubes;
    reader->AddVariable( "totalPhotons", &totalPhotons );
    reader->AddVariable( "totalCubes", &totalCubes );


    TString dir    = "dataset/weights/";
    TString prefix = "TMVARegression";
    TString weightfile = dir + prefix + "_" + methodName + ".weights.xml";
    reader->BookMVA( methodName, weightfile ); 

    TFile* tf = new TFile(inFile);
    TTree *tr = (TTree*)tf->Get("FgdLongestProjectionTree");


    tr->SetBranchAddress( "tr1", &tr1 );
    tr->SetBranchAddress( "tr2", &tr2 );
    tr->SetBranchAddress( "tr3", &tr3 );
    tr->SetBranchAddress( "ph_tr1", &ph_tr1 );
    tr->SetBranchAddress( "ph_tr2", &ph_tr2 );
    tr->SetBranchAddress( "ph_tr3", &ph_tr3 );

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
