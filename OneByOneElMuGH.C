#include "MyAnalysisElMuGH.h"
#include "Plotter.h"
#include <iostream>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include <fstream>







void MyThread(string str)
{

}

int main(int argc, char *argv[]) {

    
    
    std::string str(argv[1]);

    TString s;
    s = str;
    bool b_isData = false;
    if(s.Contains("DATA"))
    {b_isData = true;
    }
    else
    {
        b_isData = false;
    }
    std::string outFilename = str+"GH_tree.root";
    MyAnalysis* b = new  MyAnalysis(0,b_isData,outFilename);
    TChain* ch2 = new TChain("tree");
    cout <<str <<endl;
    string chstr = str+"_histoGH.root";
    ch2->Add(str.c_str());

    ch2->Process(b);
    TFile *mainf = TFile::Open(str.c_str());

    TFile* h_file =new TFile(chstr.c_str(),"RECREATE");

    h_file->cd();

    for (vector<TH1F*>::const_iterator it = b->histograms.begin(); it != b->histograms.end();++it)
    {
        (*it)->Write();


    }
    b->h_btag_eff_num->Write();
    b->h_btag_eff_den->Write();
   // b->h_btag_eff_num->Divide(b->h_btag_eff_den);
    TH1F* num_events = (TH1F*) mainf->Get("h_Nevents");
     num_events->Write("num_events");

    b->h_btag_eff_num->Write("division");

    h_file->Close();

    ch2->Delete();
    b->Delete();

	return 0;


}
