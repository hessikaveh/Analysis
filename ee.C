#include "MyAnalysis.h"
#include "Plotter.h"
#include <iostream>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <string>
#include <thread>
#include <list>
#include <TThread.h>
#include <fstream>






void MydataThread(MyAnalysis *a)
{

    TChain* ch = new TChain("tree");
    cout << "data thread" <<endl;
    ch->Add("ElEl.root");
    //        ch->Add("DoubleEGOct16ElElDATA.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1308110001.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1308430000.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1309280000.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1310000000.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1310390000.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1311110000.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1311110001.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1311500000.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1311500001.root.root");
    //        ch->Add("DoubleEGOct16ElElDATA171016_1312210000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1303520000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1304340000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1305070000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1305470000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1305470001.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1306190000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1306190001.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1306590000.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1306590001.root.root");
    //        ch->Add("SingleElectronOct16ElElDATA171016_1307300000.root.root");
    ch->Process(a);



}
void MyttThread(MyAnalysis *b)
{

    TChain* ch2 = new TChain("tree");
    cout << "tt thread ha!" <<endl;
    ch2->Add("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8TTOct16ElEl.root");

    ch2->Process(b);

}
void MydyThread(MyAnalysis *b)
{

    TChain* ch2 = new TChain("tree");
    cout << "tt thread ha!" <<endl;
    ch2->Add("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8TTOct16ElEl.root");

    ch2->Process(b);

}
void MyThread(string str)
{
    TThread::Initialize();
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
    MyAnalysis b(0,b_isData);
    TChain ch2("tree");
    cout <<str +"thread!" <<endl;
    string chstr = str+"_histo.root";
    ch2.Add(str.c_str());

    ch2.Process(&b);
    TFile h_file(chstr.c_str(),"RECREATE");

    h_file.cd();

    for (vector<TH1F*>::const_iterator it = b.histograms.begin(); it != b.histograms.end();++it)
    {
        (*it)->Write();


    }
    b.h_btag_eff_num->Write();
    b.h_btag_eff_den->Write();
    b.h_btag_eff_num->Divide(b.h_btag_eff_den);


    b.h_btag_eff_num->Write("division");

    h_file.Close();
}

int main() {
//    ROOT::EnableThreadSafety();
    unsigned int nthreads = std::thread::hardware_concurrency();
    cout << "num of threads: "<< nthreads << endl;
    float lumi = 3590;
    list<string> filenames = {"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8DY10Oct16ElElGHMC" ,
                              "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8DY10Oct16ElElMC" ,
                              "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8DY10Oct16MuMuGHMC" ,
                              "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8DY10Oct16MuMuMC"};
    vector<std::thread> t;
    vector<string> v_str;
    std::ifstream file("files.txt");
    std::string stri;
    while (std::getline(file, stri))
    {
        // Process str
        v_str.push_back(stri);
    }
    for(int i=0; i < (int) v_str.size();++i)
    {
        //        std::list<string>::iterator it = std::next(filenames.begin(), i);
        //        cout << *it << endl;
        //        std::thread t_dum(&MyThread,*it);

        std::thread t_dum(&MyThread,v_str[i]);
        cout << v_str[i]<<endl;
        t.push_back(move(t_dum));

    }
    for(vector<std::thread>::iterator it = t.begin(); it != t.end();++it)
    {
        it->join();
    }

    // Create a Proof-lite instance:
    //TProof::Open("");
    // tell the chain that we want to use PROOF
    //   chain->SetProof();
    // And this will now use all your cores!
    //   chain->Process("MySelector.C+");
    //    MyAnalysis *A = new MyAnalysis(0,true);
    //    MyAnalysis *B = new MyAnalysis(0,false);
    //    MyAnalysis *C = new MyAnalysis(0,false);

    //ch->SetProof();

    //
    //    MyttThread();
    //MydataThread();
    //    cout << "start" << endl;
    //    thread t1(&MydataThread,A);
    //    thread t2(&MyttThread,B);
    //    thread t3(&MydyThread,C);

    //    cout << "wait" <<endl;

    //    t1.join();
    //    t2.join();
    //    t3.join();
    //    cout << "done" <<endl;

    //   MyAnalysis *C = new MyAnalysis();
    //   TChain* ch3 = new TChain("events");
    //   ch3->Add("files/wjets.root");
    //   ch3->Process(C);

    //   MyAnalysis *D = new MyAnalysis();
    //   TChain* ch4 = new TChain("events");
    //   ch4->Add("files/dy.root");
    //   ch4->Process(D);

    //   MyAnalysis *E = new MyAnalysis();
    //   TChain* ch5 = new TChain("events");
    //   ch5->Add("files/ww.root");
    //   ch5->Process(E);

    //   MyAnalysis *F = new MyAnalysis();
    //   TChain* ch6 = new TChain("events");
    //   ch6->Add("files/wz.root");
    //   ch6->Process(F);

    //   MyAnalysis *G = new MyAnalysis();
    //   TChain* ch7 = new TChain("events");
    //   ch7->Add("files/zz.root");
    //   ch7->Process(G);

    //   MyAnalysis *H = new MyAnalysis();
    //   TChain* ch8 = new TChain("events");
    //   ch8->Add("files/qcd.root");
    //   ch8->Process(H);

    //   MyAnalysis *I = new MyAnalysis();
    //   TChain* ch9 = new TChain("events");
    //   ch9->Add("files/single_top.root");
    //   ch9->Process(I);


    //    TFile* h_data_file = new TFile("data_histograms.root","RECREATE");

    //    h_data_file->cd();
    //    for (vector<TH1F*>::const_iterator it = A->histograms.begin(); it != A->histograms.end();++it)
    //    {
    //        (*it)->Write();
    //    }


    //    h_data_file->Close();

    //    TFile* h_tt_file = new TFile("tt_histograms.root","RECREATE");
    //    //   TFile* h_beff_file = new TFile("beff_histograms.root","RECREATE");
    //    h_tt_file->cd();
    //    double sigTT = 831.8;
    //    cout << "b" << endl;
    //    for (vector<TH1F*>::const_iterator it = B->histograms.begin(); it != B->histograms.end();++it)
    //    {
    //        (*it)->Write();

    //        //        (*it)->Scale((lumi*sigTT)/B->h_Nevents->Integral());
    //    }
    //    B->h_btag_eff_num->Write();
    //    B->h_btag_eff_den->Write();
    //    B->h_btag_eff_num->Divide(B->h_btag_eff_den);

    //    //   h_beff_file->cd();
    //    B->h_btag_eff_num->Write("division");
    //    //       h_beff_file->Close();
    //    h_tt_file->Close();

    //    TFile* h_dy_file = new TFile("dy_histograms.root","RECREATE");
    //    h_dy_file->cd();
    //    for (vector<TH1F*>::const_iterator it = C->histograms.begin(); it != C->histograms.end();++it)
    //    {
    //        (*it)->Write();

    //        //        (*it)->Scale((lumi*sigTT)/B->h_Nevents->Integral());
    //    }
    //    h_dy_file->Close();


    //    Plotter P;
    //    P.SetData(A->histograms, std::string("Data"));
    //    P.AddBg(B->histograms, std::string("TTbar"));
    //    //	P.AddBg(C->histograms, std::string("Wjets"));
    //    //	P.AddBg(D->histograms, std::string("DY"));
    //    //	P.AddBg(E->histograms, std::string("WW"));
    //    //	P.AddBg(F->histograms, std::string("WZ"));
    //    //	P.AddBg(G->histograms, std::string("ZZ"));
    //    //	P.AddBg(H->histograms, std::string("QCD"));
    //    //	P.AddBg(I->histograms, std::string("single Top"));

    //    P.Plot(string("results.pdf"));



    //    Plotter P_MC;
    //    P_MC.AddBg(B->histograms_MC, std::string("TTbar"));
    //    //	P_MC.AddBg(C->histograms_MC, std::string("Wjets"));
    //    //	P_MC.AddBg(D->histograms_MC, std::string("DY"));
    //    //	P_MC.AddBg(E->histograms_MC, std::string("WW"));
    //    //	P_MC.AddBg(F->histograms_MC, std::string("WZ"));
    //    //	P_MC.AddBg(G->histograms_MC, std::string("ZZ"));
    //    //	P_MC.AddBg(H->histograms_MC, std::string("QCD"));
    //    //	P_MC.AddBg(I->histograms_MC, std::string("single Top"));
    //    P_MC.Plot(string("results_MC.pdf"));




}
