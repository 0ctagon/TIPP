#include <iostream>
#include <TRandom3.h>
#include "fct.cpp"

void MC()
{
    TRandom *random = new TRandom3();
    random->SetSeed(0);

    ifstream file;
    file.open("tree_struc.data");

    Particle B1;
    Particle B2;
    int nbevts = 100000;

    TFile *ftree = new TFile("MC.root","recreate");
    TTree* tree = new TTree("tree", "MC tree");
    tree->Branch("B1event", &B1, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");
    tree->Branch("B2event", &B2, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");

    RooRealVar Ep("Ep","Positron energy (GeV)",4.0);
    RooRealVar Ee("Ee","Electron energy (GeV)",7.0);
    RooRealVar angle("angle","Angle between e+e- beam (rad)",0.083);
    RooRealVar mb("mb","mass b+/b-",5.27931);
    //RooRealVar mb("mb","mass b0",5.27962);

    TLorentzVector Qp(-Ep.getVal(),0,0,Ep.getVal());
    TLorentzVector Qe(Ee.getVal()*cos(angle.getVal()),Ee.getVal()*sin(angle.getVal()),0,Ee.getVal());
    TLorentzVector Qu = Qp + Qe;

    double_t q = sqrt(Qu*Qu);
    double_t pb = sqrt(q*q-4*mb.getVal()*mb.getVal())/2;

    RooRealVar costheta("costheta","angle between B-/B in the COM",0,-1,1);
    RooRealVar phi("phi","angle between B-/B in the COM",0,0,TMath::TwoPi());

    TVector3 beta = Boost(Qu);

    for(int i=0;i<nbevts;i++)
    {
        costheta.setVal(random->Uniform(-1.0,1.0));
        phi.setVal(random->Uniform(0,TMath::TwoPi()));
        
        TLorentzVector Qb1 = TransfoLorentz(beta,costheta,phi,pb,q);
        TLorentzVector Qb2 = TransfoLorentz(beta,costheta,phi,-pb,q);

        B1.E=Qb1.E();   file >> B1.E;
        B1.P=Qb1.P();   file >> B1.P;
        B1.Px=Qb1.Px(); file >> B1.Px;
        B1.Py=Qb1.Py(); file >> B1.Py;
        B1.Pz=Qb1.Py(); file >> B1.Pz;
        B1.Costheta=costheta.getVal();file >> B1.Costheta;
        B1.Phi=phi.getVal();file >> B1.Phi;

        B2.E=Qb2.E();   file >> B2.E;
        B2.P=Qb2.P();   file >> B2.P;
        B2.Px=Qb2.Px(); file >> B2.Px;
        B2.Py=Qb2.Py(); file >> B2.Py;
        B2.Pz=Qb2.Py(); file >> B2.Pz;
        B2.Costheta=costheta.getVal();file >> B2.Costheta;
        B2.Phi=phi.getVal();file >> B2.Phi;

        tree->Fill();
        
    }
    file.close();
    tree->Write();
    ftree->Close();

}