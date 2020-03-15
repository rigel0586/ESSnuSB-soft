void tmva_1_train(TString inFile = "tmva_data.root")
{
    TFile* tf = new TFile(inFile);
    //TTree *tr = (TTree*)tf->Get("FgdTMVAData");
    TTree *tr = (TTree*)tf->Get("TotalPhotonsTree");
    

    // TClonesArray* hitPos = nullptr;
    // TClonesArray* photons = nullptr;
    // esbroot::reconstruction::superfgd::FgdTMVAEventRecord* dat = nullptr;
    // tr->SetBranchAddress("FgdTMVADataBranch", &dat);
    // tr->SetBranchAddress("FgdTMVAHitArrayBranch", &hitPos);
    // tr->SetBranchAddress("FgdTMVAPhotoArrayBranch", &photons);

    // int entries = tr->GetEntries();
    // for(int i=0; i < entries; ++i) 
    // { 
    //     cout << "               Event " << i << "           "<< endl;
    //     tr->GetEntry(i); 
    //     dat->ReadEventData();
    // }

    TMVA::Tools::Instance();

	TString outfileName( "trainTMVA.root" );
   	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
	TMVA::Factory *factory = new TMVA::Factory( "TestTMVA", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
	TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

    Double_t regWeight  = 1.0;
    dataloader->AddRegressionTree( tr, regWeight );

    // Add input variables
    dataloader->AddVariable( "totalPhotons", "Total photons", "units", 'F' );

    // Add the variable carrying the regression target
    dataloader->AddTarget( "nuEnergy" );

    factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", 
            "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut = ""; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";
    // tell the DataLoader to use all remaining events in the trees after training for testing:
    dataloader->PrepareTrainingAndTestTree( mycut,
                                         "nTrain_Regression=1000:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );

	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();


	outputFile->Close();

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> testTMVA is done!" << std::endl;

	delete factory;

	return 0;
}
