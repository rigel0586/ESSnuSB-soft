void readMCStatsExitDat()
{
    TFile* tf = new TFile("statsRoot.root");
    TTree *trExit = (TTree*)tf->Get("FgdMCExitParticlesReconstructionData");

    esbroot::reconstruction::superfgd::FgdExitData* exitDat = nullptr;
    trExit->SetBranchAddress("FgdMCExitParticlesBranch", &exitDat);


    int entries = trExit->GetEntries();
    for(int i=0; i < entries; ++i) 
    { 
        cout << "               Event " << i << "           "<< endl;
        trExit->GetEntry(i); 
        cout << "          Neutrino Pdg  " << exitDat->fnuPdg << endl;
        cout << "        Neutrino Energy " << exitDat->fnuE << endl;
        cout << "               IsWeakCC " << exitDat->fIsWeakCC << endl;
        cout << "               IsWeakNC " << exitDat->fIsWeakNC << endl;
        cout << "         IsQuasiElastic " << exitDat->fIsQuasiElastic << endl;
        const std::vector<esbroot::reconstruction::superfgd::FgdExitParticle>& exitPars =  exitDat->fVecPars;
        for(int j = 0; j < exitPars.size(); ++j)
        {
            cout << "                        Pdg " << exitPars[j].fPdg << endl;
            cout << "               Momentum.Mag " << exitPars[j].fMomentum.Mag() << endl;
            cout << "                 Momentum.X " << exitPars[j].fMomentum.X() << endl;
            cout << "                 Momentum.Y " << exitPars[j].fMomentum.Y() << endl;
            cout << "                 Momentum.Z " << exitPars[j].fMomentum.Z() << endl;
            cout << "                 Position.X " << exitPars[j].fposition.X() << endl;
            cout << "                 Position.Y " << exitPars[j].fposition.Y() << endl;
            cout << "                 Position.Z " << exitPars[j].fposition.Z() << endl;
        }
        cout << "=======================================" << endl;
    }
}


