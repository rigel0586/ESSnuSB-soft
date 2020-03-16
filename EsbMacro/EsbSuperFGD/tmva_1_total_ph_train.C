void tmva_1_total_ph_train(TString inFile = "tmva_data.root")
{
    TFile* tf = new TFile(inFile);
    TTree *tr = (TTree*)tf->Get("FgdTotalPhotonsTree");
    

    TMVA::Tools::Instance();

	TString outfileName( "trainTMVA.root" );
   	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
	TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
	TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

    // Add input variables
    dataloader->AddVariable( "totalPhotons", "Total photons", "units", 'F' );
    dataloader->AddVariable( "totalCubes", "Total cubes", "units", 'F' );

    // Add the variable carrying the regression target
    dataloader->AddTarget( "nuEnergy" );

    Double_t regWeight  = 1.0;
    dataloader->AddRegressionTree( tr, regWeight );
    

    // factory->BookMethod( dataloader,  TMVA::Types::kPDEFoam, "PDEFoam",
	// 		    "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );

    // factory->BookMethod( dataloader,  TMVA::Types::kKNN, "KNN",
    //                        "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

    // factory->BookMethod( dataloader,  TMVA::Types::kLD, "LD",
    //                        "!H:!V:VarTransform=None" );

    // factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG",
    //                        "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );

    factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut mycut = "totalPhotons>0"; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";
    // tell the DataLoader to use all remaining events in the trees after training for testing:
    dataloader->PrepareTrainingAndTestTree( mycut,
                                         "nTrain_Regression=10:nTest_Regression=5:SplitMode=Random:NormMode=NumEvents:!V" );

	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();


	outputFile->Close();

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> testTMVA is done!" << std::endl;

	delete factory;

	return 0;
}
