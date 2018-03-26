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
        return sf_nominal;
    };

    void SetJetHadFlav(int flavour)
    {
        had_flav = flavour;
    };

    int GetHadFlav()
    {
        return had_flav;
    };
    void SetJetBDis(double bDiscriminator)
    {
       b_dis = bDiscriminator;
    };
    double GetJetBDis()
    {

            return b_dis;

    };
    void SetJetBSF(double BTagSF)
    {
        b_sf = BTagSF;
    };
    double GetJetBSF()
    {
        return b_sf;
    };

    void SetJetBSFUP(double BTagSFUP)
    {
        b_sfup = BTagSFUP;
    };
    double GetJetBSFUP()
    {
        return b_sfup;
    };

    void SetJetBSFDOWN(double BTagSFDOWN)
    {
        b_sfdown = BTagSFDOWN;
    };
    double GetJetBSFDOWN()
    {
        return b_sfdown;
    };

    void SetJetPtRes(double Respt)
    {
        ptRes = Respt;
    }
    double GetJetPtRes()
    {
        return ptRes;
    }
    void SetJetPhiRes(double Resphi)
    {
        phiRes = Resphi;
    }
    double GetJetPhiRes()
    {
        return phiRes;
    }

private:

    double btag;
    bool jetid;
    TLorentzVector gen_Jet;
    double sf_nominal;
    double sf_down;
    double sf_up;
    int had_flav;
    bool b_loose,b_medium,b_tight;
    double b_dis;
    double b_sf;
    double b_sfup;
    double b_sfdown;

    double ptRes;
    double phiRes;


};

#endif /* MYJET_H_ */
