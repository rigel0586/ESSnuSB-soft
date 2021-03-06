/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers.
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables.
///
/// The methods to be used can be switched on and off by means of booleans, or
/// via the prompt command, for example:
///
///     root -l TMVARegression.C\(\"LD,MLP\"\)
///
/// (note that the backslashes are mandatory)
/// If no method given, a default set is used.
///
/// The output file "TMVAReg.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVARegression
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"


using namespace TMVA;

void tmva_train( Int_t numToTrain = 6000
                  , Int_t numToTest = 442
                  ,TString myMethodList = "" 
                  ,TString inFileName = "tmva_data.root"              // simulate_3_TMVA_data2.C and simulate_3_TMVA_data.C
                  //,TString inFileName = "tmva_graph_data.root"          // simulate_3_graph_reconstruction.C
                  //,TString inFileName = "tmva_mc_recon_data.root"     // simulate_3_MC_reconstruction.C

                  ,std::string treeName = "FgdLongestProjectionTree"      // simulate_3_TMVA_data2.C and simulate_3_TMVA_data.C 
                  //,std::string treeName = "trainTree"                     // MC muon momentum from simulate_3_graph_reconstruction.C, simulate_3_MC_reconstruction.C
                  //,std::string treeName = "fittedMomTree"                 // Genfit muon momentum from simulate_3_graph_reconstruction.C, simulate_3_MC_reconstruction.C
                  //,std::string treeName = "CalMomTree"                    // Calorimetric muon momentum from simulate_3_graph_reconstruction.C, simulate_3_MC_reconstruction.C    
   
                              )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   //     mylinux~> root -l TMVARegression.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();



   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0; // Not supported
   Use["KNN"]             = 0;
   //
   // Linear Discriminant Analysis
   Use["LD"]		        = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   //
   // Neural Network
   Use["MLP"]             = 1;
#ifdef R__HAS_TMVACPU
   Use["DNN_CPU"] = 0;
#else
   Use["DNN_CPU"] = 0;
#endif
   //
   // Support Vector Machine
   Use["SVM"]             = 0; // Not supported
   //
   // Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegression" << std::endl;

   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Create a new root output file
   TString outfileName( "TMVAReg.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );


   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //     (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //     (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   // Add input variables
   //dataloader->AddVariable( "tr1", "1st track", "units", 'F' );
   //dataloader->AddVariable( "tr2", "2nd track", "units", 'F' );
   //dataloader->AddVariable( "tr3", "3rd track", "units", 'F' );
   //dataloader->AddVariable( "ph_tr1", "1st track photons", "units", 'F' );
   //dataloader->AddVariable( "ph_tr2", "2nd track photons", "units", 'F' );
   //dataloader->AddVariable( "ph_tr3", "3rd track photons", "units", 'F' );

   //dataloader->AddVariable( "totalCubes", "Total cubes", "units", 'F' );
   //dataloader->AddVariable( "totalEdep", "Total deposited energy", "units", 'F' );

   //dataloader->AddVariable( "totalPhotons", "Total photons", "units", 'F' );
   //dataloader->AddVariable( "tr1", "1st track", "units", 'F' );
   //dataloader->AddVariable( "totalCubes", "Total cubes", "units", 'F' );

   if(treeName.compare("FgdLongestProjectionTree") !=0 )
   {
      dataloader->AddVariable( "totalPhotons", "Total photons", "units", 'F' );
      dataloader->AddVariable( "totalCubes", "Total cubes", "units", 'F' );
      dataloader->AddVariable( "muon_mom", "Muon momentum", "units", 'F' );
   }

   if(treeName.compare("FgdLongestProjectionTree") ==0 )
   {
      dataloader->AddVariable( "allEdep", "All deposited energy", "units", 'F' );
      //dataloader->AddVariable( "muon_mom", "Muon momentum", "units", 'F' );
      /*
      dataloader->AddVariable( "totalPhotons", "Total photons", "units", 'F' );
      dataloader->AddVariable( "lMean", "mean_Spectrum", "units", 'F' );
      dataloader->AddVariable( "lMeanErr", "mean_Spectrum_err", "units", 'F' );
      dataloader->AddVariable( "lStdDev", "std_dev", "units", 'F' );
      dataloader->AddVariable( "lStdDevErr", "std_dev_err", "units", 'F' );
      dataloader->AddVariable( "lEntries", "entries", "units", 'F' );
      */
      /* */
      //dataloader->AddVariable( "fPhSpecter0", "Cube photons_0", "units", 'F' ); 
      //dataloader->AddVariable( "fPhSpecter0", "Cube photons_0", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter1", "Cube photons_1", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter2", "Cube photons_2", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter3", "Cube photons_3", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter4", "Cube photons_4", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter5", "Cube photons_5", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter6", "Cube photons_6", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter7", "Cube photons_7", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter8", "Cube photons_8", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter9", "Cube photons_9", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter10", "Cube photons_10", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter11", "Cube photons_11", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter12", "Cube photons_12", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter13", "Cube photons_13", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter14", "Cube photons_14", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter15", "Cube photons_15", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter16", "Cube photons_16", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter17", "Cube photons_17", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter18", "Cube photons_18", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter19", "Cube photons_19", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter20", "Cube photons_20", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter21", "Cube photons_21", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter22", "Cube photons_22", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter23", "Cube photons_23", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter24", "Cube photons_24", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter25", "Cube photons_25", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter26", "Cube photons_26", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter27", "Cube photons_27", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter28", "Cube photons_28", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter29", "Cube photons_29", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter30", "Cube photons_30", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter31", "Cube photons_31", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter32", "Cube photons_32", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter33", "Cube photons_33", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter34", "Cube photons_34", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter35", "Cube photons_35", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter36", "Cube photons_36", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter37", "Cube photons_37", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter38", "Cube photons_38", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter39", "Cube photons_39", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter40", "Cube photons_40", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter41", "Cube photons_41", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter42", "Cube photons_42", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter43", "Cube photons_43", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter44", "Cube photons_44", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter45", "Cube photons_45", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter46", "Cube photons_46", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter47", "Cube photons_47", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter48", "Cube photons_48", "units", 'F' ); 
      dataloader->AddVariable( "fPhSpecter49", "Cube photons_49", "units", 'F' ); 
      
   }

   // Add the variable carrying the regression target
   dataloader->AddTarget( "nuEnergy" );
   //dataloader->AddTarget( "trueEdep" );
   //dataloader->AddTarget( "nuPdg" );
   //dataloader->AddTarget( "isCC" );
   //dataloader->AddTarget( "muon_mom" );

   // It is also possible to declare additional targets for multi-dimensional regression, ie:
   //     factory->AddTarget( "fvalue2" );
   // BUT: this is currently ONLY implemented for MLP

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input = TFile::Open( inFileName );
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;

   // Register the regression tree
   TTree *regTree = (TTree*)input->Get(TString(treeName));
   // global event weights per tree (see below for setting event-wise weights)
   Double_t regWeight  = 1.0;

   // You can add an arbitrary number of regression trees
   dataloader->AddRegressionTree( regTree, regWeight );

   // This would set individual event weights (the variables defined in the
   // expression need to exist in the original TTree)
   //dataloader->SetWeightExpression( "totalPhotons", "Regression" );

   // Apply additional cuts on the signal and background samples (can be different)
   //TCut mycut = "totalCubes>0 && tr1>0"; // for example: TCut mycut = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycut = "";
   // tell the DataLoader to use all remaining events in the trees after training for testing:
   TString prepareTrainStr;
   prepareTrainStr.Form("nTrain_Regression=%i:nTest_Regression=%i:SplitMode=Random:NormMode=NumEvents:!V", numToTrain, numToTest);
   dataloader->PrepareTrainingAndTestTree( mycut,prepareTrainStr);
   //
   //     dataloader->PrepareTrainingAndTestTree( mycut,
   //            "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   //     dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );

   // Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // PDE - RS method
   if (Use["PDERS"])
      factory->BookMethod( dataloader,  TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
   // And the options strings for the MinMax and RMS methods, respectively:
   //
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );

   if (Use["PDEFoam"])
       factory->BookMethod( dataloader,  TMVA::Types::kPDEFoam, "PDEFoam",
			    "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( dataloader,  TMVA::Types::kKNN, "KNN",
                           "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // Linear discriminant
   if (Use["LD"])
      factory->BookMethod( dataloader,  TMVA::Types::kLD, "LD",
                           "!H:!V:VarTransform=None" );

	// Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_MC",
                          "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=MC:SampleSize=100000:Sigma=0.1:VarTransform=D" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options) .. the formula of this example is good for parabolas
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GA",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm" );

   if (Use["FDA_MT"])
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_MT",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GAMT",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   // Neural network (MLP)
   if (Use["MLP"])
      factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BP:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

   if (Use["DNN_CPU"]) {
      /*
          TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
          TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
          TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
          TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
          TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
          TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
          TString layoutString ("Layout=TANH|100,TANH|30,LINEAR");
       */
      //TString layoutString("Layout=TANH|50,Layout=TANH|50,Layout=TANH|50,LINEAR");
      //TString layoutString("Layout=TANH|100,Layout=TANH|100,Layout=TANH|100,LINEAR");
      TString layoutString("Layout=SIGMOID|256,Layout=SIGMOID|256,Layout=SIGMOID|256,LINEAR");

      TString training0("LearningRate=1e-2,Momentum=0.5,Repetitions=1,ConvergenceSteps=20,BatchSize=100,"
                        "TestRepetitions=10,WeightDecay=0.01,Regularization=NONE,DropConfig=0.2+0.2+0.2+0.,"
                        "DropRepetitions=2");
      TString training1("LearningRate=1e-3,Momentum=0.9,Repetitions=1,ConvergenceSteps=20,BatchSize=100,"
                        "TestRepetitions=5,WeightDecay=0.01,Regularization=L2,DropConfig=0.1+0.1+0.1,DropRepetitions="
                        "1");
      TString training2("LearningRate=1e-4,Momentum=0.3,Repetitions=1,ConvergenceSteps=10,BatchSize=100,"
                        "TestRepetitions=5,WeightDecay=0.01,Regularization=NONE");

      TString trainingStrategyString("TrainingStrategy=");
      trainingStrategyString += training0 + "|" + training1 + "|" + training2;

      //       TString trainingStrategyString
      //       ("TrainingStrategy=LearningRate=1e-1,Momentum=0.3,Repetitions=3,ConvergenceSteps=20,BatchSize=30,TestRepetitions=7,WeightDecay=0.0,L1=false,DropFraction=0.0,DropRepetitions=5");

      TString nnOptions(
         "!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM:Architecture=CPU");
      //       TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
      nnOptions.Append(":");
      nnOptions.Append(layoutString);
      nnOptions.Append(":");
      nnOptions.Append(trainingStrategyString);

      factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", nnOptions); // NN
   }



   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( dataloader,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDT"])
     factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );

   if (Use["BDTG"])
     factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=10000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;

   delete factory;
   delete dataloader;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
}