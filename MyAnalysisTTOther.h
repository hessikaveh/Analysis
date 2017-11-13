//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct 21 13:11:01 2017 by ROOT version 6.08/00
// from TTree tree/tr
// found on file: DoubleEGOct15ElElDATA.root
//////////////////////////////////////////////////////////

#ifndef MyAnalysis_h
#define MyAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// Headers needed by this particular selector
#include "Math/GenVector/LorentzVector.h"

#include <vector>
#include "MyJet.h"
#include "MyLepton.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
using namespace std;

class MyAnalysis : public TSelector {
public :
    TTreeReader     fReader;  //!the tree reader
    TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

    // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> Run = {fReader, "Run"};
   TTreeReaderValue<Int_t> Event = {fReader, "Event"};
   TTreeReaderValue<Int_t> Lumi = {fReader, "Lumi"};
   TTreeReaderValue<Int_t> Bunch = {fReader, "Bunch"};
   TTreeReaderValue<Double_t> w_mc = {fReader, "w_mc"};
   TTreeReaderValue<Double_t> w_trigger_ee = {fReader, "w_trigger_ee"};
   TTreeReaderValue<Double_t> w_trigger_em = {fReader, "w_trigger_em"};
   TTreeReaderValue<Double_t> w_trigger_mm = {fReader, "w_trigger_mm"};
   TTreeReaderValue<Double_t> w_posEl = {fReader, "w_posEl"};
   TTreeReaderValue<Double_t> w_negEl = {fReader, "w_negEl"};
   TTreeReaderValue<Double_t> w_posMu = {fReader, "w_posMu"};
   TTreeReaderValue<Double_t> w_negMu = {fReader, "w_negMu"};
   TTreeReaderValue<Double_t> w_top = {fReader, "w_top"};
   TTreeReaderArray<double> w_bjets = {fReader, "w_bjets"};
   TTreeReaderValue<Float_t> rho = {fReader, "rho"};
   TTreeReaderValue<Double_t> PU_weight = {fReader, "PU_weight"};
   TTreeReaderValue<Double_t> PU_weight_secondWay = {fReader, "PU_weight_secondWay"};
   TTreeReaderArray<double> pt_Leptons = {fReader, "pt_Leptons"};
   TTreeReaderArray<double> eta_Leptons = {fReader, "eta_Leptons"};
   TTreeReaderArray<double> phi_Leptons = {fReader, "phi_Leptons"};
   TTreeReaderArray<double> e_Leptons = {fReader, "e_Leptons"};
   TTreeReaderArray<int> charge_Leptons = {fReader, "charge_Leptons"};
   TTreeReaderArray<double> pt_gen_Leptons = {fReader, "pt_gen_Leptons"};
   TTreeReaderArray<double> eta_gen_Leptons = {fReader, "eta_gen_Leptons"};
   TTreeReaderArray<double> phi_gen_Leptons = {fReader, "phi_gen_Leptons"};
   TTreeReaderArray<double> e_gen_Leptons = {fReader, "e_gen_Leptons"};
   TTreeReaderArray<int> charge_gen_Leptons = {fReader, "charge_gen_Leptons"};
   TTreeReaderArray<int> id_gen_Leptons = {fReader, "id_gen_Leptons"};
   TTreeReaderArray<double> pt_gen_daughters = {fReader, "pt_gen_daughters"};
   TTreeReaderArray<double> eta_gen_daughters = {fReader, "eta_gen_daughters"};
   TTreeReaderArray<double> phi_gen_daughters = {fReader, "phi_gen_daughters"};
   TTreeReaderArray<double> e_gen_daughters = {fReader, "e_gen_daughters"};
   TTreeReaderArray<int> id_gen_daughters = {fReader, "id_gen_daughters"};
   TTreeReaderArray<double> pt_gen_mothers = {fReader, "pt_gen_mothers"};
   TTreeReaderArray<double> eta_gen_mothers = {fReader, "eta_gen_mothers"};
   TTreeReaderArray<double> phi_gen_mothers = {fReader, "phi_gen_mothers"};
   TTreeReaderArray<double> e_gen_mothers = {fReader, "e_gen_mothers"};
   TTreeReaderArray<int> id_gen_mothers = {fReader, "id_gen_mothers"};
   TTreeReaderArray<double> pt_Jets = {fReader, "pt_Jets"};
   TTreeReaderArray<double> eta_Jets = {fReader, "eta_Jets"};
   TTreeReaderArray<double> phi_Jets = {fReader, "phi_Jets"};
   TTreeReaderArray<double> e_Jets = {fReader, "e_Jets"};
   TTreeReaderArray<int> hadflav_Jets = {fReader, "hadflav_Jets"};
   TTreeReaderArray<double> ptresolution_Jets = {fReader, "ptresolution_Jets"};
   TTreeReaderArray<double> phiresolution_Jets = {fReader, "phiresolution_Jets"};
   TTreeReaderArray<double> sf_nominal_Jets = {fReader, "sf_nominal_Jets"};
   TTreeReaderArray<double> sf_up_Jets = {fReader, "sf_up_Jets"};
   TTreeReaderArray<double> sf_down_Jets = {fReader, "sf_down_Jets"};
   TTreeReaderValue<vector<bool>> Loose_Bdiscriminator = {fReader, "Loose_Bdiscriminator"};
   TTreeReaderValue<vector<bool>> Medium_Bdiscriminator = {fReader, "Medium_Bdiscriminator"};
   TTreeReaderValue<vector<bool>> Tight_Bdiscriminator = {fReader, "Tight_Bdiscriminator"};
   TTreeReaderArray<double> pt_gen_Jets = {fReader, "pt_gen_Jets"};
   TTreeReaderArray<double> eta_gen_Jets = {fReader, "eta_gen_Jets"};
   TTreeReaderArray<double> phi_gen_Jets = {fReader, "phi_gen_Jets"};
   TTreeReaderArray<double> e_gen_Jets = {fReader, "e_gen_Jets"};
   TTreeReaderArray<double> hadFlavour_gen_Jets = {fReader, "hadFlavour_gen_Jets"};
   TTreeReaderArray<double> partonFlavour_gen_Jets = {fReader, "partonFlavour_gen_Jets"};
   TTreeReaderArray<double> pt_bJets = {fReader, "pt_bJets"};
   TTreeReaderArray<double> eta_bJets = {fReader, "eta_bJets"};
   TTreeReaderArray<double> phi_bJets = {fReader, "phi_bJets"};
   TTreeReaderArray<double> e_bJets = {fReader, "e_bJets"};
   TTreeReaderArray<double> pt_mets = {fReader, "pt_mets"};
   TTreeReaderArray<double> px_mets = {fReader, "px_mets"};
   TTreeReaderArray<double> py_mets = {fReader, "py_mets"};
   TTreeReaderValue<Int_t> num_PU_vertices = {fReader, "num_PU_vertices"};
   TTreeReaderValue<Int_t> PU_BunchCrossing = {fReader, "PU_BunchCrossing"};
   TTreeReaderValue<Int_t> num_PU_gen_vertices = {fReader, "num_PU_gen_vertices"};
   TTreeReaderValue<Int_t> num_PV = {fReader, "num_PV"};


    MyAnalysis(TTree * /*tree*/ =0,bool isData=0){b_isData = isData;}
    virtual ~MyAnalysis() { }
    virtual Int_t   Version() const { return 2; }
    virtual void    Begin(TTree *tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Notify();
    virtual Bool_t  Process(Long64_t entry);
    virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void    SetOption(const char *option) { fOption = option; }
    virtual void    SetObject(TObject *obj) { fObject = obj; }
    virtual void    SetInputList(TList *input) { fInput = input; }
    virtual TList  *GetOutputList() const { return fOutput; }
    virtual void    SlaveTerminate();
    virtual void    Terminate();
    //Functions and Variables for inputing the data into classes
    void BuildEvents();
    void bEffCalc(double b_eta, double b_pt, vector<MyJet> jets);
    bool applySF(bool isBTagged, double Btag_SF, double Btag_eff);

    //Defining Histograms
    int TotalEvents =0;
    double weight = 1;
    vector<TH1F*> histograms;
    vector<TH1F*> histograms_MC;
    TH1F* h_Nevents;
    TH1F* h_Nevents_ALS;
    TH1F* h_Nevents_AJS;
    TH1F* h_Nevents_AMS;
    TH1F* h_Nevents_ABS;

    TH1F* h_num_PV;
    TH1F* h_num_PV_weighted;
    TH1F* h_m_dilepton_ALS_NoVeto;
    TH1F* h_num_jets_ALS_NoVeto;
    TH1F* h_pt_lepton_ALS;
    TH1F* h_MET_ALS;
    TH1F* h_num_in_noMET;
    TH1F* h_num_in_MET;
    TH1F* h_num_out_noMET;
    TH1F* h_num_out_MET;
    TH1F* h_Jets_pt;
    TH1F* h_Jets_eta;
    TH1F* h_Jets_num;
    TH1F* h_gen_Jets_pt;
    TH1F* h_gen_Jets_eta;
    TH1F* h_Leptons_pt;
    TH2F* h_btag_eff;
    TH2F* h_btag_eff_num;
    TH2F* h_btag_eff_den;
    TH1F* h_btag_ef;
    TH1F* h_N_btagged_b;
    TH1F* h_N_btagged_total;
    vector<MyJet> MyJets;
    vector<MyJet> BJets;
    vector<MyLepton> MyLeptons;
    double met_pt;
    //ClassDef(MyAnalysis,0);
    BTagCalibrationReader reader;
    string btagSf="CSVv2_Moriond17_B_H.csv";
    TFile* f_btagEff;
    bool b_isData;
};

#endif

#ifdef MyAnalysis_cxx
void MyAnalysis::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    fReader.SetTree(tree);
}

Bool_t MyAnalysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}


#endif // #ifdef MyAnalysis_cxx
