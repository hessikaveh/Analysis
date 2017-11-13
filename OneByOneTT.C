#include "MyAnalysisTT.h"
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
    MyAnalysis* b = new  MyAnalysis(0,b_isData);
    TChain* ch2 = new TChain("tree");
    cout <<str <<endl;
    string chstr = str+"_histoTT.root";
    ch2->Add(str.c_str());

    ch2->Process(b);
    TFile* h_file =new TFile(chstr.c_str(),"RECREATE");

    h_file->cd();

    for (vector<TH1F*>::const_iterator it = b->histograms.begin(); it != b->histograms.end();++it)
    {
        (*it)->Write();


    }
    b->h_btag_eff_num->Write();
    b->h_btag_eff_den->Write();
    b->h_btag_eff_num->Divide(b->h_btag_eff_den);


    b->h_btag_eff_num->Write("division");

    h_file->Close();



	return 0;


}
