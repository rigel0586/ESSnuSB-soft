void tmva_2_long_tr_app(TString inFile = "tmva_data.root", TString methodName = "BDTG")
{
    TMVA::Tools::Instance();
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 

    Float_t x_tr0, y_tr0, z_tr0, x_ph_tr0, y_ph_tr0, z_ph_tr0;
    reader->AddVariable( "x_tr0", &x_tr0 );
    reader->AddVariable( "y_tr0", &y_tr0 );
    reader->AddVariable( "z_tr0", &z_tr0 );
    reader->AddVariable( "x_ph_tr0", &x_ph_tr0 );
    reader->AddVariable( "y_ph_tr0", &y_ph_tr0 );
    reader->AddVariable( "z_ph_tr0", &z_ph_tr0 );


    Float_t x_tr1, y_tr1, z_tr1, x_ph_tr1, y_ph_tr1, z_ph_tr1;
    reader->AddVariable( "x_tr1", &x_tr1 );
    reader->AddVariable( "y_tr1", &y_tr1 );
    reader->AddVariable( "z_tr1", &z_tr1 );
    reader->AddVariable( "x_ph_tr1", &x_ph_tr1 );
    reader->AddVariable( "y_ph_tr1", &y_ph_tr1 );
    reader->AddVariable( "z_ph_tr1", &z_ph_tr1 );

    Float_t x_tr2, y_tr2, z_tr2, x_ph_tr2, y_ph_tr2, z_ph_tr2;
    reader->AddVariable( "x_tr2", &x_tr2 );
    reader->AddVariable( "y_tr2", &y_tr2 );
    reader->AddVariable( "z_tr2", &z_tr2 );
    reader->AddVariable( "x_ph_tr2", &x_ph_tr2 );
    reader->AddVariable( "y_ph_tr2", &y_ph_tr2 );
    reader->AddVariable( "z_ph_tr2", &z_ph_tr2 );

    Float_t totalPhotons, totalCubes;
    reader->AddVariable( "totalPhotons", &totalPhotons );
    reader->AddVariable( "totalCubes", &totalCubes );


    TString dir    = "dataset/weights/";
    TString prefix = "TMVARegression";
    TString weightfile = dir + prefix + "_" + methodName + ".weights.xml";
    reader->BookMVA( methodName, weightfile ); 

    TFile* tf = new TFile(inFile);
    TTree *tr = (TTree*)tf->Get("FgdLongestProjectionTree");


    tr->SetBranchAddress( "x_tr0", &x_tr0 );
    tr->SetBranchAddress( "y_tr0", &y_tr0 );
    tr->SetBranchAddress( "z_tr0", &z_tr0 );
    tr->SetBranchAddress( "x_ph_tr0", &x_ph_tr0 );
    tr->SetBranchAddress( "y_ph_tr0", &y_ph_tr0 );
    tr->SetBranchAddress( "z_ph_tr0", &z_ph_tr0 );

    tr->SetBranchAddress( "x_tr1", &x_tr1 );
    tr->SetBranchAddress( "y_tr1", &y_tr1 );
    tr->SetBranchAddress( "z_tr1", &z_tr1 );
    tr->SetBranchAddress( "x_ph_tr1", &x_ph_tr1 );
    tr->SetBranchAddress( "y_ph_tr1", &y_ph_tr1 );
    tr->SetBranchAddress( "z_ph_tr1", &z_ph_tr1 );

    tr->SetBranchAddress( "x_tr2", &x_tr2 );
    tr->SetBranchAddress( "y_tr2", &y_tr2 );
    tr->SetBranchAddress( "z_tr2", &z_tr2 );
    tr->SetBranchAddress( "x_ph_tr2", &x_ph_tr2 );
    tr->SetBranchAddress( "y_ph_tr2", &y_ph_tr2 );
    tr->SetBranchAddress( "z_ph_tr2", &z_ph_tr2 );

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
