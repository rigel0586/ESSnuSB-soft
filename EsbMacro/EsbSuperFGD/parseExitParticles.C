void parseExitParticles(TString inFile = "/*** TODO PATH ***/statsRoot.root"            // Path to statsRoot.root generated from simulate_3_MC_lepton_stats.C
                        , TString outFile = "/*** TODO PATH ***/formattedData.txt"      // output file where to write the nuance formatted data
                        , Bool_t printData = true)                                      // true - to print the info to the console
{
    TFile* tf = new TFile(inFile);
    TTree *tr = (TTree*)tf->Get("FgdMCExitParticlesReconstructionData");

    esbroot::reconstruction::superfgd::FgdExitData* dat = nullptr;
    tr->SetBranchAddress("FgdMCExitParticlesBranch", &dat);

    std::ofstream outputFile(outFile, std::ios::trunc);

    if(outputFile.is_open())
	{
        

        int entries = tr->GetEntries();
        for(int i=0; i < entries; ++i) 
        { 
            outputFile << "$ begin"  << endl;
            outputFile << "$ nuance 0"  << endl;

            tr->GetEntry(i); 
            
            if(printData)
            {
                cout << "               Event " << i << "           "<< endl;
                cout << "           Neutrino pdg " << dat->fnuPdg << endl;
                cout << "        Neutrino Energy " << dat->fnuE << endl;
                cout << "               IsWeakCC " << dat->fIsWeakCC<< endl;
                cout << "               IsWeakNC " << dat->fIsWeakNC << endl;
                cout << "         IsQuasiElastic " << dat->fIsQuasiElastic << endl;
            }
            

            std::vector<esbroot::reconstruction::superfgd::FgdExitParticle>& pars = dat->fVecPars;
            for(int j = 0; j < pars.size(); ++j)
            {
                esbroot::reconstruction::superfgd::FgdExitParticle& p = pars[j];
                if(printData)
                {
                    cout << "                        Pdg " << p.fPdg << endl;
                    cout << "                   Momentum " << p.fMomentum.X() << " " << p.fMomentum.Y() << " " << p.fMomentum.Z() << endl;
                    cout << "                   Position " << p.fposition.X() << " " << p.fposition.Y() << " " << p.fposition.Z() << endl;
                }
                // Print Vertex
                outputFile << "$ info 0 " << (j/10) << " " << (j%10) << endl;
                outputFile << "$ vertex " << p.fposition.X() << " " << p.fposition.Y() << " " << p.fposition.Z() << " 0 " << endl;
                outputFile << "$ track " << p.fPdg << " " << p.fMomentum.Mag() << " "
                           << p.fMomentum.X() << " " << p.fMomentum.Y() << " " << p.fMomentum.Z() << " 0 " << endl;
            }

            outputFile << "$ end"  << endl;
            
            if(printData) {cout << "=======================================" << endl;}
        }

        outputFile << "$ stop"  << endl;
    }
    outputFile.close();
}
