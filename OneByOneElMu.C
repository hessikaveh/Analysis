#include "MyAnalysisElMu.h"
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
    std::string str2(argv[2]);
    TString s;
    s = str;
    TString s2;
    s2 = str2;
    bool b_isData = false;
    if(s.Contains("DATA"))
    {b_isData = true;
    }
    else
    {
        b_isData = false;
    }
    bool b_isGH = false;
    bool b_isEE = false;
    bool b_isEM = false;
    bool b_isMM = false;

    if(s2.Contains("ee")) b_isEE = true;
    if(s2.Contains("mm")) b_isMM = true;
    if(s2.Contains("em")) b_isEM = true;
    if(s2.Contains("GH")) b_isGH = true;
    std::string outFilename;
     string chstr;
    if(b_isGH)
    {
       outFilename  = "out_dir/"+str+"_treeGH.root";
       chstr = "out_dir/"+str+"_histoGH.root";
    }
    else
    {
         outFilename = "out_dir/"+str+"_tree.root";
        chstr = "out_dir/"+str+"_histo.root";
    }
    MyAnalysis* b = new  MyAnalysis(0,b_isData,b_isGH,b_isEE,b_isMM,b_isEM,outFilename);
    TChain* ch2 = new TChain("tree");
    cout << outFilename << chstr <<endl;

    ch2->Add(str.c_str());

    ch2->Process(b);
    TFile *mainf = TFile::Open(str.c_str());
    TFile* h_file =new TFile(chstr.c_str(),"RECREATE");

    h_file->cd();

    for (vector<TH1F*>::const_iterator it = b->histograms.begin(); it != b->histograms.end();++it)
    {
        (*it)->Write();


    }

    TH1F* num_events = (TH1F*) mainf->Get("h_Nevents");
    num_events->Write("num_events");
    b->h_btag_eff_num->Write();
    b->h_btag_eff_den->Write();
    // b->h_btag_eff_num->Divide(b->h_btag_eff_den);


    b->h_btag_eff_num->Write("division");

    h_file->Close();
    ch2->Delete();
    b->Delete();



    return 0;


}
