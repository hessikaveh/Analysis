#define MyAnalysis_cxx
// The class definition in MyAnalysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("MyAnalysis.C")
// root> T->Process("MyAnalysis.C","some options")
// root> T->Process("MyAnalysis.C+")
//


#include "MyAnalysisElMuGHTT.h"
#include <TH2.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TGraph.h>


void MyAnalysis::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    //Histogram Definitions

    h_Nevents = new TH1F("h_Nevents","Number of events",1,0,2);
    h_Nevents->Sumw2();
    histograms.push_back(h_Nevents);
    histograms_MC.push_back(h_Nevents);

    h_Nevents_ALS = new TH1F("h_Nevents_ALS","Number of Events ALS",1,0,2);
    h_Nevents_ALS->Sumw2();
    histograms.push_back(h_Nevents_ALS);
    histograms_MC.push_back(h_Nevents_ALS);

    h_Nevents_AJS = new TH1F("h_Nevents_AJS","Number of Events AJS",1,0,2);
    h_Nevents_AJS->Sumw2();
    histograms.push_back(h_Nevents_AJS);
    histograms_MC.push_back(h_Nevents_AJS);

    h_Nevents_AMS = new TH1F("h_Nevents_AMS","Number of Events AMS",1,0,2);
    h_Nevents_AMS->Sumw2();
    histograms.push_back(h_Nevents_AMS);
    histograms_MC.push_back(h_Nevents_AMS);

    h_Nevents_ABS = new TH1F("h_Nevents_ABS","Number of Events ABS",1,0,2);
    h_Nevents_ABS->Sumw2();
    histograms.push_back(h_Nevents_ABS);
    histograms_MC.push_back(h_Nevents_ABS);

    h_num_PV = new TH1F("h_num_PV","Number of Pileup vertices",50,0,50);
    h_num_PV->Sumw2();
    histograms.push_back(h_num_PV);
    histograms_MC.push_back(h_num_PV);

    h_num_PV_weighted = new TH1F("h_num_PV_weighted","Number of Pileup vertices",50,0,50);
    h_num_PV_weighted->Sumw2();
    histograms.push_back(h_num_PV_weighted);
    histograms_MC.push_back(h_num_PV_weighted);

    h_m_dilepton_ALS_NoVeto = new TH1F("h_m_dilepton_ALS_NoVeto","Invariant mass of dilepton",50,0,300);
    h_m_dilepton_ALS_NoVeto->Sumw2();
    histograms.push_back(h_m_dilepton_ALS_NoVeto);
    histograms_MC.push_back(h_m_dilepton_ALS_NoVeto);

    h_num_jets_ALS_NoVeto = new TH1F("h_num_jets_ALS_NoVeto","Number of Jets",9,0,9);
    h_m_dilepton_ALS_NoVeto->Sumw2();
    histograms.push_back(h_num_jets_ALS_NoVeto);
    histograms_MC.push_back(h_num_jets_ALS_NoVeto);

    h_num_in_MET = new TH1F("h_num_in_MET","",1,0,2);
    h_num_in_MET->Sumw2();
    histograms.push_back(h_num_in_MET);
    histograms_MC.push_back(h_num_in_MET);

    h_num_in_noMET = new TH1F("h_num_in_noMET","",1,0,2);
    h_num_in_noMET->Sumw2();
    histograms.push_back(h_num_in_noMET);
    histograms_MC.push_back(h_num_in_noMET);

    h_num_out_MET = new TH1F("h_num_out_MET","",1,0,2);
    h_num_out_MET->Sumw2();
    histograms.push_back(h_num_out_MET);
    histograms_MC.push_back(h_num_out_MET);

    h_num_out_noMET = new TH1F("h_num_out_noMET","",1,0,2);
    h_num_out_noMET->Sumw2();
    histograms.push_back(h_num_out_noMET);
    histograms_MC.push_back(h_num_out_noMET);

    h_pt_lepton_ALS = new TH1F("h_pt_lepton_ALS","pt of leptons",38,25,400);
    h_pt_lepton_ALS->Sumw2();
    histograms.push_back(h_pt_lepton_ALS);
    histograms_MC.push_back(h_pt_lepton_ALS);

    h_MET_ALS = new TH1F("h_MET_ALS","MET",30,0,300);
    h_MET_ALS->Sumw2();
    histograms.push_back(h_MET_ALS);
    histograms_MC.push_back(h_MET_ALS);

    h_MET_sum_ALS = new TH1F("h_MET_sum_ALS","MET_sum",30,0,300);
    h_MET_sum_ALS->Sumw2();
    histograms.push_back(h_MET_sum_ALS);
    histograms_MC.push_back(h_MET_sum_ALS);

    h_Jets_pt = new TH1F("h_Jets_pt","Pt of Jets",35,0,300);
    h_Jets_pt->Sumw2();
    histograms.push_back(h_Jets_pt);
    histograms_MC.push_back(h_Jets_pt);

    h_gen_Jets_pt = new TH1F("h_gen_Jets_pt","Pt of gen Jets",35,0,300);
    h_gen_Jets_pt->Sumw2();
    histograms.push_back(h_gen_Jets_pt);
    histograms_MC.push_back(h_gen_Jets_pt);

    h_Jets_eta = new TH1F("h_Jets_eta","Eta of Jets",35,-3,3);
    h_Jets_eta->Sumw2();
    histograms.push_back(h_Jets_eta);
    histograms_MC.push_back(h_Jets_eta);

    h_gen_Jets_eta = new TH1F("h_gen_Jets_eta","Eta of gen_Jets",35,-3,3);
    h_gen_Jets_eta->Sumw2();
    histograms.push_back(h_gen_Jets_eta);
    histograms_MC.push_back(h_gen_Jets_eta);

    h_Jets_num = new TH1F("h_Jets_num","num of Jets",8,0,8);
    h_Jets_num->Sumw2();
    histograms.push_back(h_Jets_num);
    histograms_MC.push_back(h_Jets_num);

    h_Leptons_pt = new TH1F("h_Leptons_pt","Pt of Leptons",35,0,300);
    h_Leptons_pt->Sumw2();
    histograms.push_back(h_Leptons_pt);
    histograms_MC.push_back(h_Leptons_pt);

    h_btag_ef = new TH1F("h_btag_ef","efficiency",100,0,5);
    h_btag_ef->Sumw2();
    histograms.push_back(h_btag_ef);
    histograms_MC.push_back(h_btag_ef);
    h_btag_eff = new TH2F("h_btag_eff","b-tagging eff",100,0,1000,5,0,2.5);
    h_btag_eff->Sumw2();
    h_btag_eff_num = new TH2F("h_btag_eff_num","b-tagging eff num",50,30,500,5,0,2.5);
    h_btag_eff_num->Sumw2();
    h_btag_eff_den = new TH2F("h_btag_eff_den","b-tagging eff den",50,30,500,5,0,2.5);
    h_btag_eff_den->Sumw2();

    BTagCalibration calib("csvv2", btagSf);
    reader = BTagCalibrationReader(BTagEntry::OP_LOOSE,  // operating point
                                   "central",             // central sys type
    {"up", "down"}
                                   );      // other sys types

    reader.load(calib,                // calibration instance
                BTagEntry::FLAV_B,    // btag flavour
                "comb");              // measurement type
    f_btagEff = new TFile("beff_histograms.root");

    if(!b_isData){

        f_muonIDSF = new TFile(muonIDSF.c_str());
        f_muonISOSF = new TFile(muonISOSF.c_str());

        f_muonTkSF = new TFile(muonTkSF.c_str());

        f_mm_sf = new TFile(mm_sf.c_str());
        f_egammaSF = new TFile(egammaSF.c_str());
        f_egammaTkSF = new TFile(egammaTkSF.c_str());
        f_ee_sf = new TFile(ee_sf.c_str());

        f_em_sf = new TFile(em_sf.c_str());

    }
    TString option = GetOption();
}

void MyAnalysis::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

}

Bool_t MyAnalysis::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // When processing keyed objects with PROOF, the object is already loaded
    // and is available via the fObject pointer.
    //
    // This function should contain the \"body\" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.
    ++TotalEvents;
    fReader.SetEntry(entry);
    weight = 1.;
    if (TotalEvents % 1000 == 0)
        cout << "Next event -----> " << TotalEvents << endl;
    //            if(TotalEvents > 100000) return kTRUE;
    BuildEvents();
    met_pt =  *(px_mets.begin()) + *(py_mets.begin());
    met_pt *=-1;
    met_sum = *(sum_mets.begin());
    met_sum *= -1;
    // cout << b_isData << endl;
    if(!b_isData)weight *= (*w_mc);
    if(MyLeptons.size() < 2) return kTRUE;
    MyLepton lep1 = MyLeptons.at(0);
    MyLepton lep2 = MyLeptons.at(1);
    

    h_Nevents->Fill(1);


    if(MyLeptons.at(0).Pt() > MyLeptons.at(1).Pt())
    {
        if(MyLeptons.at(0).Pt() < 25) return kTRUE;
    }
    else
    {
        if(MyLeptons.at(1).Pt() < 25) return kTRUE;
    }

    double m_dilepton =-999;
    m_dilepton = (MyLeptons.at(0)+MyLeptons.at(1)).M();
    if(m_dilepton <= 20 ) return kTRUE;
    if(!(fabs(lep1.GetGenLepDaughterId()) == 24 && fabs(lep1.GetGenLepMotherId()) == 6 )) return kTRUE;
     if(!(fabs(lep2.GetGenLepDaughterId()) == 24 && fabs(lep2.GetGenLepMotherId()) == 6 )) return kTRUE;
    
    if(m_dilepton < 20 ) cout << "cut why?!"<<endl;
    h_num_PV->Fill(*num_PV,weight);
    if(!b_isData) weight *= (*PU_weight);
    h_num_PV_weighted->Fill(*num_PV,weight);
    if(!b_isData){

        TH2F* h2D_muonIDSF = (TH2F*) f_muonIDSF->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
        TH2F* h2D_muonISOSF = (TH2F*) f_muonISOSF->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
        TGraph* g_muonTkSF = (TGraph*) f_muonTkSF->Get("ratio_eff_aeta_dr030e030_corr");
        TH2F* h2D_mm_sf = (TH2F*) f_mm_sf->Get("scalefactor_eta2d_with_syst");

        TH2F* h2D_egammaSF = (TH2F*) f_egammaSF->Get("EGamma_SF2D");
        TH2F* h2D_egammaTkSF = (TH2F*) f_egammaTkSF->Get("EGamma_SF2D");

        for(vector<MyLepton>::iterator it = MyLeptons.begin(); it != MyLeptons.end();++it)
        {
            if(it->GetLepType() == "muon"){
                Int_t binx = h2D_muonIDSF->GetXaxis()->FindBin(fabs(it->Eta()));
                Int_t biny = h2D_muonIDSF->GetYaxis()->FindBin(it->Pt());
                Int_t binx2 = h2D_muonISOSF->GetXaxis()->FindBin(fabs(it->Eta()));
                Int_t biny2 = h2D_muonISOSF->GetYaxis()->FindBin(fabs(it->Pt()));
                (it->Pt() > 119) ? biny = h2D_muonIDSF->GetYaxis()->FindBin(118):biny = h2D_muonIDSF->GetYaxis()->FindBin(it->Pt());
                (it->Pt() > 119) ? biny2 = h2D_muonISOSF->GetYaxis()->FindBin(118):biny2 = h2D_muonISOSF->GetYaxis()->FindBin(it->Pt());
                Double_t TkSF = g_muonTkSF->Eval(fabs(it->Eta()));
                //cout <<"Mupt: "<< posMu.pt() << " eta: " << fabs(posMu.eta()) <<" SF: "<< h2D_muonIDSF->GetBinContent(binx,biny)<<"& "<< h2D_muonISOSF->GetBinContent(binx2,biny2) <<" & "<< TkSF <<endl;
                double SF_posMu = h2D_muonIDSF->GetBinContent(binx,biny)*h2D_muonISOSF->GetBinContent(binx2,biny2)*TkSF;
                weight *= SF_posMu;
            }
            if(it->GetLepType() == "electron")
            {
                Int_t binx = h2D_egammaSF->GetXaxis()->FindBin(it->GetLepSCeta());
                Int_t biny = h2D_egammaSF->GetYaxis()->FindBin(it->Pt());
                Int_t binx2 = h2D_egammaTkSF->GetXaxis()->FindBin(it->GetLepSCeta());
                Int_t biny2 = 1;
                double SF_posEl = h2D_egammaSF->GetBinContent(binx,biny)*h2D_egammaTkSF->GetBinContent(binx2,biny2);

                weight *= SF_posEl;
            }

        }
        TH2F* h2D_em_sf = (TH2F*) f_em_sf->Get("scalefactor_eta2d_with_syst");
        Int_t binx = h2D_em_sf->GetXaxis()->FindBin(fabs(lep1.Eta()));
        Int_t biny = h2D_em_sf->GetYaxis()->FindBin(fabs(lep2.Eta()));
        double d_em_sf = h2D_em_sf->GetBinContent(binx,biny);
        weight *= d_em_sf;



    }
    //if(!b_isData) weight *= (*w_trigger_ee) * (*w_trigger_em) *(*w_trigger_mm)*(*w_posEl)*(*w_negEl)*(*w_posMu)*(*w_negMu);
    h_m_dilepton_ALS_NoVeto->Fill(m_dilepton,weight);
    h_num_jets_ALS_NoVeto->Fill(MyJets.size(),weight);
    if(m_dilepton >= 76 && m_dilepton <= 106)
    {
        if(MyJets.size() >= 2 )
        {
            h_num_in_noMET->Fill(1,weight);
            if(met_pt >= 40) h_num_in_MET->Fill(1,weight);
        }
    }
    if(m_dilepton >= 76 && m_dilepton <= 106) return kTRUE;
    h_Nevents_ALS->Fill(1,weight);
    if(MyJets.size() >= 2 )
    {
        h_num_out_noMET->Fill(1,weight);
        if(met_sum > 40) h_num_out_MET->Fill(1,weight);
    }

    for(vector<MyLepton>::iterator it = MyLeptons.begin(); it != MyLeptons.end();++it)
    {
        h_pt_lepton_ALS->Fill(it->Pt(),weight);
        h_Leptons_pt->Fill(it->Pt(),weight);
    }
    h_MET_sum_ALS->Fill(met_sum,weight);
    h_MET_ALS->Fill(met_pt,weight);
    if(MyJets.size()<2) return kTRUE;
    double w_jet = 1;

    double bweight=1.;
    TH2F* h2D_btagEff = (TH2F*) f_btagEff->Get("division");
    if(!b_isData){
        for(vector<MyJet>::iterator it = MyJets.begin(); it != MyJets.end();++it)
        {
            Int_t binx = h2D_btagEff->GetXaxis()->FindBin(it->Pt());
            Int_t biny = h2D_btagEff->GetYaxis()->FindBin(fabs(it->Eta()));
            bEffCalc(it->Eta(),it->Pt(),MyJets);
            if(applySF(it->GetJetBDis(0),it->GetJetBSF(),h2D_btagEff->GetBinContent(binx,biny))) bweight*= it->GetJetBSF();
        }
    }
    for(vector<MyJet>::iterator it = MyJets.begin(); it != MyJets.end();++it)
    {



        //       h_btag_eff->Fill(it->Pt(),fabs(it->Eta()),N_btagged_b/N_btagged_total);
        w_jet*= it->GetJetSF(1);
        h_Jets_pt->Fill(it->Pt(),weight);
        if(it->GetGenJet().Pt()>0) h_gen_Jets_pt->Fill(it->GetGenJet().Pt(),weight);
        h_Jets_eta->Fill(it->Eta(),weight);
        if(it->GetGenJet().Pt()>0) h_gen_Jets_eta->Fill(it->GetGenJet().Eta(),weight);
        if(it->GetJetBDis(0))
        {
            BJets.push_back(*it);
        }


    }
    //    h_btag_eff->SetBinContent(h_btag_eff->FindBin(5),h_btag_eff->FindBin(fabs(.2)),N_btagged_b/N_btagged_total);

    h_Nevents_AJS->Fill(1,weight);

    h_Jets_num->Fill(MyJets.size(),weight);

    //cout << met_pt << endl;
    if(met_sum <= 40) return kTRUE;
    h_Nevents_AMS->Fill(1,weight);

    if(BJets.size() < 1) return kTRUE;
    if(!b_isData) weight*=bweight;
    h_Nevents_ABS->Fill(1,weight);
    return kTRUE;
}

void MyAnalysis::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void MyAnalysis::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
    //    top_pt->Draw("hist");
    //    top_pt->Write("hist.pdf");

}

void MyAnalysis::BuildEvents()
{
    MyJets.clear();
    MyLeptons.clear();
    vector<MyLepton> MyLeptonsD;
    MyLeptonsD.clear();
    vector<MyLepton> MyLeptonsE;
    MyLeptonsE.clear();
    vector<MyLepton> MyLeptonsM;
    MyLeptonsE.clear();
    for(int i=0; i < (int) pt_Jets.GetSize(); ++i)
    {
        MyJet jet = MyJet(pt_Jets[i],eta_Jets[i],phi_Jets[i],e_Jets[i]);
        jet.SetJetSF(sf_down_Jets[i],sf_nominal_Jets[i],sf_up_Jets[i]);
        jet.SetJetBDis(Loose_Bdiscriminator->at(i),Medium_Bdiscriminator->at(i),Tight_Bdiscriminator->at(i));
        jet.SetJetHadFlav(hadflav_Jets[i]);
        double jet_scalefactor    = reader.eval_auto_bounds(
                    "central",
                    BTagEntry::FLAV_B,
                    eta_Jets[i],
                    pt_Jets[i]
                    );
        jet.SetJetBSF(jet_scalefactor);
        if(!(pt_gen_Jets.IsEmpty())) jet.SetGenJet(pt_gen_Jets[i],eta_gen_Jets[i],phi_gen_Jets[i],e_gen_Jets[i]);
        MyJets.push_back(jet);
    }

    for(int j=0; j< (int) pt_Leptons.GetSize();++j)
    {
        if(type_Leptons[j] != "electron") continue;


        MyLepton lep = MyLepton(pt_Leptons[j],eta_Leptons[j],phi_Leptons[j],e_Leptons[j]);
lep.SetLepType("electron");
        if(!(pt_gen_Leptons.IsEmpty())) lep.SetGenLep(pt_gen_Leptons[j],eta_gen_Leptons[j],phi_gen_Leptons[j],e_gen_Leptons[j]);
        lep.SetGenLepDaughterId(id_gen_daughters[j]);
        lep.SetGenLepMotherId(id_gen_mothers[j]);
        MyLeptonsD.push_back(lep);
    }
    for(int j=0; j< (int) pt_Leptons.GetSize();++j)
    {
        if(type_Leptons[j] != "muon") continue;


        MyLepton lep = MyLepton(pt_Leptons[j],eta_Leptons[j],phi_Leptons[j],e_Leptons[j]);
        lep.SetLepType("muon");
        if(!(pt_gen_Leptons.IsEmpty())) lep.SetGenLep(pt_gen_Leptons[j],eta_gen_Leptons[j],phi_gen_Leptons[j],e_gen_Leptons[j]);
        lep.SetGenLepDaughterId(id_gen_daughters[j]);
        lep.SetGenLepMotherId(id_gen_mothers[j]);
        MyLeptonsM.push_back(lep);
    }
    for(int j=0; j< (int) MyLeptonsD.size();++j)
    {
        MyLeptonsD.at(j).SetLepSCeta(SuperClusterEta_Leptons[j]);
        if(fabs(SuperClusterEta_Leptons[j]) <= 1.479){
            if(!(d0_Leptons[j] < 0.05) && !(dz_Leptons[j] < 0.10) ) continue;
        }

        if(fabs(SuperClusterEta_Leptons[j]) > 1.479 ){
            if(!(d0_Leptons[j] < 0.10) && !(dz_Leptons[j] < 0.20) ) continue;
        }
        MyLeptonsE.push_back(MyLeptonsD.at(j));
    }
    if(!(MyLeptonsE.size() > 0 && MyLeptonsM.size() >0)) return;

    for(int i=0; i< (int) (MyLeptonsE.size()+MyLeptonsM.size());++i)
    {
        if((MyLeptonsE.size()+MyLeptonsM.size()) < 2) continue;
        if(MyLeptonsE.at(0).Pt() > MyLeptonsM.at(0).Pt())
        {
            MyLeptons.push_back(MyLeptonsE.at(0));
            MyLeptons.push_back(MyLeptonsM.at(0));
        }
        if(MyLeptonsE.at(0).Pt() < MyLeptonsM.at(0).Pt())
        {

            MyLeptons.push_back(MyLeptonsM.at(0));
            MyLeptons.push_back(MyLeptonsE.at(0));
        }
    }

}

void MyAnalysis::bEffCalc(double b_eta, double b_pt,vector<MyJet> jets)
{


    for(vector<MyJet>::iterator it = jets.begin(); it != jets.end();++it)
    {

        if(it->GetHadFlav() == 5 && it->GetJetBDis(0))  h_btag_eff_num->Fill(b_pt,fabs(b_eta));
        if(it->GetJetBDis(0)) h_btag_eff_den->Fill(b_pt,fabs(b_eta));

    }





}

bool MyAnalysis::applySF(bool isBTagged, double Btag_SF, double Btag_eff)
{
    bool newBTag = isBTagged;
    TRandom* rand_ = new TRandom;
    if (Btag_SF == 1) return newBTag; //no correction needed

    //throw die

    float coin = rand_->Uniform(1.);

    if(Btag_SF > 1){  // use this if SF>1

        if( !isBTagged ) {

            //fraction of jets that need to be upgraded
            float mistagPercent = (1.0 - Btag_SF) / (1.0 - (1.0/Btag_eff) );

            //upgrade to tagged
            if( coin < mistagPercent ) {newBTag = true;}
        }

    }else{  // use this if SF<1

        //downgrade tagged to untagged
        if( isBTagged && coin > Btag_SF ) {newBTag = false;}

    }

    return newBTag;

}
