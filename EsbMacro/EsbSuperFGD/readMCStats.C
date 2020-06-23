void readMCStats()
{
    TFile* tf = new TFile("statsRoot.root");
    TTree *tr = (TTree*)tf->Get("FgdMCLeptonStatsReconstructionData");

    esbroot::reconstruction::superfgd::FgdMCEventRecord* dat = nullptr;
    tr->SetBranchAddress("FgdMCLeptonStatsBranch", &dat);

    std::ofstream outputFile("../../EsbMacro/tests/dedx.dat", std::ios::trunc);
    if(outputFile.is_open())
    {
        outputFile << "#dE/dx" << "  " <<  "muon momentum " << endl;
    }

    double interval = 0.1;
    int bins = 50 / interval;
    std::map<int ,std::vector<double>> vals;
    for(size_t i = 0; i < 35; ++i)
        vals[i] = std::vector<double>();

    int entries = tr->GetEntries();
    for(int i=0; i < entries; ++i) 
    { 
        cout << "               Event " << i << "           "<< endl;
        tr->GetEntry(i); 
        dat->ReadEventData();
        const std::vector<std::pair<Int_t, TVector3>>& primPars =  dat->GetPrimaryParticles();
        cout << "          Muon Momentum " << dat->GetMuonMom().Mag() << endl;
        cout << "        Neutrino Energy " << dat->GetNuE() << endl;
        cout << "    IsPrimaryLeptonMuon " << dat->IsPrimaryLeptonMuon() << endl;
        cout << "               IsWeakCC " << dat->IsWeakCC() << endl;
        cout << "               IsWeakNC " << dat->IsWeakNC() << endl;
        cout << "         IsQuasiElastic " << dat->IsQuasiElastic() << endl;
        cout << "                   DeDx " << dat->GetDeDx() << endl;
        for(int j = 0; j < primPars.size(); ++j)
        {
            cout << "                        Pdg " << primPars[j].first << endl;
        }
        cout << "=======================================" << endl;

        // if(dat->IsPrimaryLeptonMuon() && outputFile.is_open()
        //     && dat->GetDeDx()!=0)
        // {
        //     outputFile <<  dat->GetDeDx()  << "  " <<  dat->GetMuonMom().Mag()<< endl;
        // }

        // if(dat->IsPrimaryLeptonMuon() && outputFile.is_open()
        //     && dat->GetEdepDeDx()!=0)
        // {
        //     outputFile <<  dat->GetMuonTrackLengthOrigin() * dat->GetEdepDeDx()  << "  " <<  dat->GetMuonMom().Mag()<< endl;
        // }

        int key = dat->GetDeDx() / interval;
        vals[key].emplace_back(dat->GetMuonMom().Mag());
    }

    std::map<int ,std::vector<double>>::iterator it = vals.begin();

    while(it != vals.end())
    {
        std::vector<double>& vec = it->second;
        double avg(0);
        for(size_t i=0; i < vec.size(); ++i)
        {
            avg+=vec[i];
        }
        avg= avg/vec.size();

        double stddev(0);
        for(size_t i=0; i < vec.size(); ++i)
        {
            stddev = (vec[i] - avg) * (vec[i] - avg);
        }

        int del = vec.size() == 1 ? 1 : vec.size() - 1;
        stddev = stddev/del;
        stddev = std::sqrt(stddev);

        outputFile <<  it->first*interval  << "  " <<  avg << "  " <<  stddev << endl;
        ++it;
    }



    outputFile.close();
}



