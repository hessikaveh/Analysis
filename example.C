#include "MyAnalysis.h"
#include "Plotter.h"
#include <iostream>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <string>
//#include <TProof.h>
int main() {
   
   float lumi = 3590;
   // Create a Proof-lite instance:
//   TProof::Open("");
   // tell the chain that we want to use PROOF
//   chain->SetProof();
   // And this will now use all your cores!
//   chain->Process("MySelector.C+");
   MyAnalysis *A = new MyAnalysis(0,true);

   TChain* ch = new TChain("tree");
   ch->Add("files/data.root");
   ch->Process(A);

   MyAnalysis *B = new MyAnalysis(0,false);
   TChain* ch2 = new TChain("tree");
   ch2->Add("files/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8TTOct15ElElMC.root");
   ch2->Process(B);
   
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
   TFile* h_data_file = new TFile("data_histograms.root","RECREATE");

   h_data_file->cd();
   for (vector<TH1F*>::const_iterator it = A->histograms.begin(); it != A->histograms.end();++it)
   {
       (*it)->Write();
   }
   A->h_btag_eff->Write();
   A->h_btag_eff_num->Write();
   A->h_btag_eff_den->Write();
   A->h_btag_eff_num->Divide(A->h_btag_eff_den);
   A->h_btag_eff_num->Write("division");
   h_data_file->Close();

   TFile* h_tt_file = new TFile("tt_histograms.root","RECREATE");
//   TFile* h_beff_file = new TFile("beff_histograms.root","RECREATE");
   h_tt_file->cd();
   double sigTT = 831.8;
   for (vector<TH1F*>::const_iterator it = B->histograms.begin(); it != B->histograms.end();++it)
   {
       (*it)->Write();
       (*it)->Scale((lumi*sigTT)/B->h_Nevents->Integral());
   }
   B->h_btag_eff_num->Write();
   B->h_btag_eff_den->Write();
   B->h_btag_eff_num->Divide(B->h_btag_eff_den);
   h_tt_file->Close();
//   h_beff_file->cd();
   B->h_btag_eff_num->Write("division");
//   h_beff_file->Close();

	Plotter P;
	P.SetData(A->histograms, std::string("Data"));
	P.AddBg(B->histograms, std::string("TTbar"));
//	P.AddBg(C->histograms, std::string("Wjets"));
//	P.AddBg(D->histograms, std::string("DY"));
//	P.AddBg(E->histograms, std::string("WW"));
//	P.AddBg(F->histograms, std::string("WZ"));
//	P.AddBg(G->histograms, std::string("ZZ"));
//	P.AddBg(H->histograms, std::string("QCD"));
//	P.AddBg(I->histograms, std::string("single Top"));
   
	P.Plot(string("results.pdf"));


   
	Plotter P_MC;
	P_MC.AddBg(B->histograms_MC, std::string("TTbar"));
//	P_MC.AddBg(C->histograms_MC, std::string("Wjets"));
//	P_MC.AddBg(D->histograms_MC, std::string("DY"));
//	P_MC.AddBg(E->histograms_MC, std::string("WW"));
//	P_MC.AddBg(F->histograms_MC, std::string("WZ"));
//	P_MC.AddBg(G->histograms_MC, std::string("ZZ"));
//	P_MC.AddBg(H->histograms_MC, std::string("QCD"));
//	P_MC.AddBg(I->histograms_MC, std::string("single Top"));
   P_MC.Plot(string("results_MC.pdf"));

}
