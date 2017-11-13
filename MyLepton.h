/*
 * MyLepton.h
 *
 *  Created on: Feb 1, 2012
 *      Author: csander
 */

#ifndef MYLepton_H_
#define MYLepton_H_

#include <TLorentzVector.h>

class MyLepton: public TLorentzVector {

public:

    MyLepton();
    MyLepton(double pt, double eta, double phi, double e) {
        SetPtEtaPhiE(pt, eta, phi, e);
    }
    ;
    virtual ~MyLepton();

    void SetIsolation(double x) {
        isolation = x;
    }
    ;

    void SetCharge(int q) {
        charge = q;
    }
    ;

    const double GetIsolation() {
        return isolation;
    }
    ;
    const bool IsIsolated() {
        return (isolation < 1.);
    }
    ;
    const int GetCharge() {
        return charge;
    }
    ;
    void SetGenLep(double gen_pt,double gen_eta,double gen_phi,double gen_e)
    {
        gen_lep.SetPtEtaPhiE(gen_pt,gen_eta,gen_phi,gen_e);
    }

    TLorentzVector GetGenLep()
    {
        return gen_lep;
    }
    void SetGenLepMotherId(int id)
    {
        id_gen_mother = id;
    }
    void SetGenLepDaughterId(int id)
    {
        id_gen_daughter = id;
    }
    int GetGenLepMotherId()
    {
        return id_gen_mother;
    }
    int GetGenLepDaughterId()
    {
        return id_gen_daughter;
    }
    void SetLepSCeta(double SCeta)
    {
        eta_SC = SCeta;
    }
    void SetLepd0(double lepd0)
    {
        d0 = lepd0;
    }
    void SetLepdZ(double lepdZ)
    {
        dZ = lepdZ;
    }
    double GetLepSCeta()
    {
        return eta_SC;
    }
    double GetLepd0()
    {
        return d0;
    }
    double GetLepdZ()
    {
        return dZ;
    }
    void SetLepType(std::string type)
    {
        lep_type = type;
    }
    std::string GetLepType()
    {
        return lep_type;
    }
private:

    double isolation;
    int charge;
    TLorentzVector gen_lep;
    int id_gen_daughter;
    int id_gen_mother;
    double eta_SC;
    double d0;
    double dZ;
    std::string lep_type;

};

#endif /* MYLepton_H_ */
