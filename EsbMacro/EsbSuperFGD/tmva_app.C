/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides a simple example on how to use the trained regression MVAs
/// within an analysis module
///
///  - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
///  - Package   : TMVA
///  - Exectuable: TMVARegressionApplication
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;

void tmva_app( TString myMethodList = ""
               ,TString inFileName = "tmva_data.root"              // simulate_3_TMVA_data2.C and simulate_3_TMVA_data.C
               //,TString inFileName = "tmva_graph_data.root"          // simulate_3_graph_reconstruction.C
               //,TString inFileName = "tmva_mc_recon_data.root"     // simulate_3_MC_reconstruction.C

               ,std::string treeName = "FgdLongestProjectionTree"      // simulate_3_TMVA_data2.C and simulate_3_TMVA_data.C 
               //,std::string treeName = "trainTree"                       // MC muon momentum from simulate_3_graph_reconstruction.C, simulate_3_MC_reconstruction.C
               //,std::string treeName = "fittedMomTree"                 // Genfit muon momentum from simulate_3_graph_reconstruction.C, simulate_3_MC_reconstruction.C
               //,std::string treeName = "CalMomTree"                    // Calorimetric muon momentum from simulate_3_graph_reconstruction.C, simulate_3_MC_reconstruction.C    
   
               ) 
{
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0; // Not supported
   Use["KNN"]             = 0;
   // 
   // --- Linear Discriminant Analysis
   Use["LD"]		        = 0;
   // 
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   // 
   // --- Neural Network
   Use["MLP"] = 1;
#ifdef R__HAS_TMVACPU
   Use["DNN_CPU"] = 0;
#else
   Use["DNN_CPU"] = 0;
#endif

   // 
   // --- Support Vector Machine 
   Use["SVM"]             = 0; // Not supported
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegressionApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

    Float_t tr1, tr2, tr3, ph_tr1, ph_tr2, ph_tr3;
    //reader->AddVariable( "tr1", &tr1 );
    //reader->AddVariable( "tr2", &tr2 );
    //reader->AddVariable( "tr3", &tr3 );
    //reader->AddVariable( "ph_tr1", &ph_tr1 );
    //reader->AddVariable( "ph_tr2", &ph_tr2 );
    //reader->AddVariable( "ph_tr3", &ph_tr3 );

    Float_t totalPhotons, totalCubes, totalEdep, allEdep, muon_mom;
    //reader->AddVariable( "totalPhotons", &totalPhotons );
    //reader->AddVariable( "totalCubes", &totalCubes );
    //reader->AddVariable( "totalEdep", &totalEdep );

    Float_t ph_sp_0, ph_sp_1, ph_sp_2, ph_sp_3, ph_sp_4, ph_sp_5, ph_sp_6, ph_sp_7, ph_sp_8, ph_sp_9, 
            ph_sp_10, ph_sp_11, ph_sp_12, ph_sp_13, ph_sp_14, ph_sp_15, ph_sp_16, ph_sp_17, ph_sp_18, ph_sp_19, 
            ph_sp_20, ph_sp_21, ph_sp_22, ph_sp_23, ph_sp_24, ph_sp_25, ph_sp_26, ph_sp_27, ph_sp_28, ph_sp_29, 
            ph_sp_30, ph_sp_31, ph_sp_32, ph_sp_33, ph_sp_34, ph_sp_35, ph_sp_36, ph_sp_37, ph_sp_38, ph_sp_39, 
            ph_sp_40, ph_sp_41, ph_sp_42, ph_sp_43, ph_sp_44, ph_sp_45, ph_sp_46, ph_sp_47, ph_sp_48, ph_sp_49;

   Float_t lMean, lMeanErr, lStdDev, lStdDevErr, lEntries;

    //reader->AddVariable( "totalPhotons", &totalPhotons );
    //reader->AddVariable( "tr1", &tr1 );
    //reader->AddVariable( "totalCubes", &totalCubes );
    if(treeName.compare("FgdLongestProjectionTree") !=0 )
    {
      reader->AddVariable( "totalPhotons", &totalPhotons );
      reader->AddVariable( "totalCubes", &totalCubes );
      reader->AddVariable( "muon_mom", &muon_mom );
    }

    if(treeName.compare("FgdLongestProjectionTree") ==0 )
    {
         reader->AddVariable( "allEdep", &allEdep );
         //reader->AddVariable( "muon_mom", &muon_mom );
         /* 
         reader->AddVariable( "totalPhotons", &totalPhotons );
         reader->AddVariable( "lMean", &lMean );
         reader->AddVariable( "lMeanErr", &lMeanErr );
         reader->AddVariable( "lStdDev", &lStdDev );
         reader->AddVariable( "lStdDevErr", &lStdDevErr );
         reader->AddVariable( "lEntries", &lEntries );
         */

         //reader->AddVariable( "fPhSpecter0", &ph_sp_0 ); 
         reader->AddVariable( "fPhSpecter1", &ph_sp_1 ); 
         reader->AddVariable( "fPhSpecter2", &ph_sp_2 ); 
         reader->AddVariable( "fPhSpecter3", &ph_sp_3 ); 
         reader->AddVariable( "fPhSpecter4", &ph_sp_4 ); 
         reader->AddVariable( "fPhSpecter5", &ph_sp_5 ); 
         reader->AddVariable( "fPhSpecter6", &ph_sp_6 ); 
         reader->AddVariable( "fPhSpecter7", &ph_sp_7 ); 
         reader->AddVariable( "fPhSpecter8", &ph_sp_8 ); 
         reader->AddVariable( "fPhSpecter9", &ph_sp_9 ); 
         reader->AddVariable( "fPhSpecter10", &ph_sp_10 ); 
         reader->AddVariable( "fPhSpecter11", &ph_sp_11 ); 
         reader->AddVariable( "fPhSpecter12", &ph_sp_12 ); 
         reader->AddVariable( "fPhSpecter13", &ph_sp_13 ); 
         reader->AddVariable( "fPhSpecter14", &ph_sp_14 ); 
         reader->AddVariable( "fPhSpecter15", &ph_sp_15 ); 
         reader->AddVariable( "fPhSpecter16", &ph_sp_16 ); 
         reader->AddVariable( "fPhSpecter17", &ph_sp_17 ); 
         reader->AddVariable( "fPhSpecter18", &ph_sp_18 ); 
         reader->AddVariable( "fPhSpecter19", &ph_sp_19 ); 
         reader->AddVariable( "fPhSpecter20", &ph_sp_20 ); 
         reader->AddVariable( "fPhSpecter21", &ph_sp_21 ); 
         reader->AddVariable( "fPhSpecter22", &ph_sp_22 ); 
         reader->AddVariable( "fPhSpecter23", &ph_sp_23 ); 
         reader->AddVariable( "fPhSpecter24", &ph_sp_24 ); 
         reader->AddVariable( "fPhSpecter25", &ph_sp_25 ); 
         reader->AddVariable( "fPhSpecter26", &ph_sp_26 ); 
         reader->AddVariable( "fPhSpecter27", &ph_sp_27 ); 
         reader->AddVariable( "fPhSpecter28", &ph_sp_28 ); 
         reader->AddVariable( "fPhSpecter29", &ph_sp_29 ); 
         reader->AddVariable( "fPhSpecter30", &ph_sp_30 ); 
         reader->AddVariable( "fPhSpecter31", &ph_sp_31 ); 
         reader->AddVariable( "fPhSpecter32", &ph_sp_32 ); 
         reader->AddVariable( "fPhSpecter33", &ph_sp_33 ); 
         reader->AddVariable( "fPhSpecter34", &ph_sp_34 ); 
         reader->AddVariable( "fPhSpecter35", &ph_sp_35 ); 
         reader->AddVariable( "fPhSpecter36", &ph_sp_36 ); 
         reader->AddVariable( "fPhSpecter37", &ph_sp_37 ); 
         reader->AddVariable( "fPhSpecter38", &ph_sp_38 ); 
         reader->AddVariable( "fPhSpecter39", &ph_sp_39 ); 
         reader->AddVariable( "fPhSpecter40", &ph_sp_40 ); 
         reader->AddVariable( "fPhSpecter41", &ph_sp_41 ); 
         reader->AddVariable( "fPhSpecter42", &ph_sp_42 ); 
         reader->AddVariable( "fPhSpecter43", &ph_sp_43 ); 
         reader->AddVariable( "fPhSpecter44", &ph_sp_44 ); 
         reader->AddVariable( "fPhSpecter45", &ph_sp_45 ); 
         reader->AddVariable( "fPhSpecter46", &ph_sp_46 ); 
         reader->AddVariable( "fPhSpecter47", &ph_sp_47 ); 
         reader->AddVariable( "fPhSpecter48", &ph_sp_48 ); 
         reader->AddVariable( "fPhSpecter49", &ph_sp_49 ); 
       
    }

   // --- Book the MVA methods

   TString dir    = "dataset/weights/";
   TString prefix = "TMVARegression";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
         TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   TH1* hists[100];
   Int_t nhists = -1;
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      TH1* h = new TH1F( it->first.c_str(), TString(it->first) + " method", 100, 0, 2 );
      if (it->second) hists[++nhists] = h;
   }
   nhists++;
   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input(0);
   if (!gSystem->AccessPathName( inFileName )) {
      input = TFile::Open( inFileName ); // check if file in local directory exists
   }
   else {
      std::cout << "ERROR: Input file is not accessible" << std::endl;
      exit(1);
   }
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegressionApp        : Using input file: " << input->GetName() << std::endl;

   // --- Event loop

   // Prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   TTree *theTree = (TTree*)input->Get(TString(treeName));    

   std::cout << "--- Select signal sample" << std::endl;
   //theTree->SetBranchAddress( "tr1", &tr1 );
   //theTree->SetBranchAddress( "tr2", &tr2 );
   //theTree->SetBranchAddress( "tr3", &tr3 );
   //theTree->SetBranchAddress( "ph_tr1", &ph_tr1 );
   //theTree->SetBranchAddress( "ph_tr2", &ph_tr2 );
   //theTree->SetBranchAddress( "ph_tr3", &ph_tr3 );

   if(treeName.compare("FgdLongestProjectionTree") !=0 )
   {
      theTree->SetBranchAddress( "totalPhotons", &totalPhotons );
      theTree->SetBranchAddress( "totalCubes", &totalCubes );
      //theTree->SetBranchAddress( "totalEdep", &totalEdep );
      theTree->SetBranchAddress( "muon_mom", &muon_mom );
   }

   if(treeName.compare("FgdLongestProjectionTree") ==0 )
   {
      theTree->SetBranchAddress( "allEdep", &allEdep );
      //theTree->SetBranchAddress( "muon_mom", &muon_mom );
      /*
      theTree->SetBranchAddress( "totalPhotons", &totalPhotons );
      theTree->SetBranchAddress( "lMean", &lMean );
      theTree->SetBranchAddress( "lMeanErr", &lMeanErr );
      theTree->SetBranchAddress( "lStdDev", &lStdDev );
      theTree->SetBranchAddress( "lStdDevErr", &lStdDevErr );
      theTree->SetBranchAddress( "lEntries", &lEntries );
      */

      //theTree->SetBranchAddress( "fPhSpecter0", &ph_sp_0 ); 
      theTree->SetBranchAddress( "fPhSpecter1", &ph_sp_1 ); 
      theTree->SetBranchAddress( "fPhSpecter2", &ph_sp_2 ); 
      theTree->SetBranchAddress( "fPhSpecter3", &ph_sp_3 ); 
      theTree->SetBranchAddress( "fPhSpecter4", &ph_sp_4 ); 
      theTree->SetBranchAddress( "fPhSpecter5", &ph_sp_5 ); 
      theTree->SetBranchAddress( "fPhSpecter6", &ph_sp_6 ); 
      theTree->SetBranchAddress( "fPhSpecter7", &ph_sp_7 ); 
      theTree->SetBranchAddress( "fPhSpecter8", &ph_sp_8 ); 
      theTree->SetBranchAddress( "fPhSpecter9", &ph_sp_9 ); 
      theTree->SetBranchAddress( "fPhSpecter10", &ph_sp_10 ); 
      theTree->SetBranchAddress( "fPhSpecter11", &ph_sp_11 ); 
      theTree->SetBranchAddress( "fPhSpecter12", &ph_sp_12 ); 
      theTree->SetBranchAddress( "fPhSpecter13", &ph_sp_13 ); 
      theTree->SetBranchAddress( "fPhSpecter14", &ph_sp_14 ); 
      theTree->SetBranchAddress( "fPhSpecter15", &ph_sp_15 ); 
      theTree->SetBranchAddress( "fPhSpecter16", &ph_sp_16 ); 
      theTree->SetBranchAddress( "fPhSpecter17", &ph_sp_17 ); 
      theTree->SetBranchAddress( "fPhSpecter18", &ph_sp_18 ); 
      theTree->SetBranchAddress( "fPhSpecter19", &ph_sp_19 ); 
      theTree->SetBranchAddress( "fPhSpecter20", &ph_sp_20 ); 
      theTree->SetBranchAddress( "fPhSpecter21", &ph_sp_21 ); 
      theTree->SetBranchAddress( "fPhSpecter22", &ph_sp_22 ); 
      theTree->SetBranchAddress( "fPhSpecter23", &ph_sp_23 ); 
      theTree->SetBranchAddress( "fPhSpecter24", &ph_sp_24 ); 
      theTree->SetBranchAddress( "fPhSpecter25", &ph_sp_25 ); 
      theTree->SetBranchAddress( "fPhSpecter26", &ph_sp_26 ); 
      theTree->SetBranchAddress( "fPhSpecter27", &ph_sp_27 ); 
      theTree->SetBranchAddress( "fPhSpecter28", &ph_sp_28 ); 
      theTree->SetBranchAddress( "fPhSpecter29", &ph_sp_29 ); 
      theTree->SetBranchAddress( "fPhSpecter30", &ph_sp_30 ); 
      theTree->SetBranchAddress( "fPhSpecter31", &ph_sp_31 ); 
      theTree->SetBranchAddress( "fPhSpecter32", &ph_sp_32 ); 
      theTree->SetBranchAddress( "fPhSpecter33", &ph_sp_33 ); 
      theTree->SetBranchAddress( "fPhSpecter34", &ph_sp_34 ); 
      theTree->SetBranchAddress( "fPhSpecter35", &ph_sp_35 ); 
      theTree->SetBranchAddress( "fPhSpecter36", &ph_sp_36 ); 
      theTree->SetBranchAddress( "fPhSpecter37", &ph_sp_37 ); 
      theTree->SetBranchAddress( "fPhSpecter38", &ph_sp_38 ); 
      theTree->SetBranchAddress( "fPhSpecter39", &ph_sp_39 ); 
      theTree->SetBranchAddress( "fPhSpecter40", &ph_sp_40 ); 
      theTree->SetBranchAddress( "fPhSpecter41", &ph_sp_41 ); 
      theTree->SetBranchAddress( "fPhSpecter42", &ph_sp_42 ); 
      theTree->SetBranchAddress( "fPhSpecter43", &ph_sp_43 ); 
      theTree->SetBranchAddress( "fPhSpecter44", &ph_sp_44 ); 
      theTree->SetBranchAddress( "fPhSpecter45", &ph_sp_45 ); 
      theTree->SetBranchAddress( "fPhSpecter46", &ph_sp_46 ); 
      theTree->SetBranchAddress( "fPhSpecter47", &ph_sp_47 ); 
      theTree->SetBranchAddress( "fPhSpecter48", &ph_sp_48 ); 
      theTree->SetBranchAddress( "fPhSpecter49", &ph_sp_49 ); 
   }

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++)
   {
        if (ievt%1000 == 0) {
            std::cout << "--- ... Processing event: " << ievt << std::endl;
        }

        theTree->GetEntry(ievt);

        // Retrieve the MVA target values (regression outputs) and fill into histograms
        // NOTE: EvaluateRegression(..) returns a vector for multi-target regression

        for (Int_t ih=0; ih<nhists; ih++) {
            TString title = hists[ih]->GetTitle();
            Float_t val = (reader->EvaluateRegression( title ))[0];
            hists[ih]->Fill( val);    
        }
   }
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   //canvas->SaveAs("nu_E");

   // --- Write histograms
   TFile *target  = new TFile( "TMVARegApp.root","RECREATE" );
   for (Int_t ih=0; ih<nhists; ih++)
   {
       hists[ih]->Write();
   } 
   target->Close();

   std::cout << "--- Created root file: \"" << target->GetName() 
             << "\" containing the MVA output histograms" << std::endl;
  
   delete reader;
    
   std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}