/*
 * MyJet.h
 *
 *  Created on: Feb 1, 2012
 *      Author: csander
 */

#ifndef MYJET_H_
#define MYJET_H_

#include <TLorentzVector.h>
#include <TRandom3.h>

class MyJet: public TLorentzVector {

public:

    MyJet();
    MyJet(double pt, double eta, double phi, double e) {
        SetPtEtaPhiE(pt, eta, phi, e);
    }

    virtual ~MyJet();

    void SetGenJet(double gen_pt,double gen_eta,double gen_phi,double gen_e)
    {
        gen_Jet.SetPtEtaPhiE(gen_pt,gen_eta,gen_phi,gen_e);
    };

    TLorentzVector GetGenJet()
    {
        return gen_Jet;
    };
    void SetJetSF(double down,double nominal,double up)
    {
        sf_down = down;
        sf_nominal = nominal;
        sf_up = up;

    };
    double GetJetSF(int i)
    {
        if(i == 1){
            return sf_nominal;
        }
        else if(i == 0)
        {
            return sf_down;
        }
        else if(i == 2)
        {
            return sf_up;
        }
        return 1.;
    };

    void SetJetHadFlav(int flavour)
    {
        had_flav = flavour;
    };

    int GetHadFlav()
    {
        return had_flav;
    };
    void SetJetBDis(bool loose, bool medium, bool tight)
    {
        b_loose = loose;
        b_medium = medium;
        b_tight = tight;
    };
    bool GetJetBDis(int i)
    {
        if(i==0)
        {
            return b_loose;
        }
        else if(i == 1)
        {
            return b_medium;
        }
        else if(i==2 )
        {
            return b_tight;
        }
        return 1.;
    };
    void SetJetBSF(double BTagSF)
    {
        b_sf = BTagSF;
    };
    double GetJetBSF()
    {
        return b_sf;
    };

private:

    double btag;
    bool jetid;
    TLorentzVector gen_Jet;
    double sf_nominal;
    double sf_down;
    double sf_up;
    int had_flav;
    bool b_loose,b_medium,b_tight;
    double b_sf;


};

#endif /* MYJET_H_ */
