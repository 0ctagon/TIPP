#include <iostream>
#include <TRandom3.h>
#include "fct.cpp"

void MC()
{
    //Set seed
    TRandom *random = new TRandom3();
    random->SetSeed(0);

    ifstream file;
    file.open("tree_struc.data");

    //Define particles for our simulation
    Particle B1;
    Particle B2;
    Particle K;
    Particle Nu1;
    Particle Nu2;
    Particle NuNu;

    int nbevts = 100000;

    //Define the tree, each branch is a particle
    TFile *ftree = new TFile("MC.root","recreate");
    TTree* tree = new TTree("tree", "MC tree");
    tree->Branch("B1event", &B1, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");
    tree->Branch("B2event", &B2, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");
    tree->Branch("Kevent", &K, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");
    tree->Branch("Nu1event", &Nu1, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");
    tree->Branch("Nu2event", &Nu2, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");
    tree->Branch("NuNuevent", &NuNu, "E/D:P/D:Px/D:Py/D:Pz/D:Costheta/D:Phi/D");

    //All the constant in our simulation
    RooRealVar Ep("Ep","Positron energy (GeV)",4.0);
    RooRealVar Ee("Ee","Electron energy (GeV)",7.0);
    RooRealVar angle("angle","Angle between e+e- beam (rad)",0.083);
    RooRealVar mb("mb","mass b+/b-",5.27931);
    RooRealVar mk("mk","mass k+",0.493677);
    //RooRealVar mb("mb","mass b0",5.27962);
    //RooRealVar mk("mk","mass k0",0.497611);

    //Define the 4-vectors for e-/e+ & Upsilon 4S
    TLorentzVector Qp(-Ep.getVal(),0,0,Ep.getVal());
    TLorentzVector Qe(Ee.getVal()*cos(angle.getVal()),-Ee.getVal()*sin(angle.getVal()),0,Ee.getVal());
    TLorentzVector Qu = Qp + Qe;

    //Calculate q and the impulsion of the B mesons
    double_t q = sqrt(Qu*Qu);
    double_t pb = sqrt(q*q-4*mb.getVal()*mb.getVal())/2;

    //Definition of all the random variables
    RooRealVar costheta("costheta","angle between B-/B in the COM",0,-1,1);
    RooRealVar phi("phi","angle between B-/B in the COM",0,0,TMath::TwoPi());

    RooRealVar costheta1("costheta1","angle between K/2nu in the COM of B",0,-1,1);
    RooRealVar phi1("phi","angle between K/2nu in the COM of B",0,0,TMath::TwoPi());

    RooRealVar costheta2("costheta2","angle between 2nu in the COM of B",0,-1,1);
    RooRealVar phi2("phi2","angle between 2nu in the COM of B",0,0,TMath::TwoPi());

    RooRealVar m12("m12","mass 2nu system",0,0,(mb.getVal()-mk.getVal()));

    //Calculate the boost for the B mesons
    TVector3 betaU = Boost(Qu);

    for(int i=0;i<nbevts;i++)
    {
        //Get values for all the random variable
        costheta.setVal(random->Uniform(costheta.getMin(),costheta.getMax()));
        phi.setVal(random->Uniform(phi.getMin(),phi.getMax()));

        costheta1.setVal(random->Uniform(costheta1.getMin(),costheta1.getMax()));
        phi1.setVal(random->Uniform(phi1.getMin(),phi1.getMax()));

        costheta2.setVal(random->Uniform(costheta2.getMin(),costheta2.getMax()));
        phi2.setVal(random->Uniform(phi2.getMin(),phi2.getMax()));
        
        m12.setVal(random->Uniform(m12.getMin(),m12.getMax()));
        
        //Calculate the 4-vector for B & B_
        TLorentzVector Qb1 = TransfoLorentz(betaU,costheta,phi,pb,q/2);
        TLorentzVector Qb2 = TransfoLorentz(betaU,costheta,phi,-pb,q/2);

        //Calculate the 4-vector for K & 2nu
        TVector3 betaB1 = Boost(Qb1);

        double_t pK = sqrt((pow(mb.getVal(),2)-pow(m12.getVal()+mk.getVal(),2))*(pow(mb.getVal(),2)-pow(m12.getVal()-mk.getVal(),2)))/(2*mb.getVal());
        double_t EK = sqrt(mk.getVal()*mk.getVal()+pK*pK);
        
        double_t pnunu = -pK;
        double_t Enunu = (q/2)-EK;

        TLorentzVector Qk = TransfoLorentz(betaB1, costheta1, phi1, pK, EK);
        TLorentzVector Qnunu = TransfoLorentz(betaB1, costheta1, phi1, pnunu, Enunu);

        //Calculate the 4-vector for nu1 & nu2
        TVector3 betanunu = Boost(Qnunu);

        double_t qnunu = sqrt(Qnunu*Qnunu);

        TLorentzVector Qnu1 = TransfoLorentz(betanunu, costheta2, phi2, qnunu/2, qnunu/2);
        TLorentzVector Qnu2 = TransfoLorentz(betanunu, costheta2, phi2, -qnunu/2, qnunu/2);

        //Building the tree
        B1.E=Qb1.E();   file >> B1.E;
        B1.P=Qb1.P();   file >> B1.P;
        B1.Px=Qb1.Px(); file >> B1.Px;
        B1.Py=Qb1.Py(); file >> B1.Py;
        B1.Pz=Qb1.Pz(); file >> B1.Pz;
        B1.Costheta=costheta.getVal();file >> B1.Costheta;
        B1.Phi=phi.getVal();file >> B1.Phi;

        B2.E=Qb2.E();   file >> B2.E;
        B2.P=Qb2.P();   file >> B2.P;
        B2.Px=Qb2.Px(); file >> B2.Px;
        B2.Py=Qb2.Py(); file >> B2.Py;
        B2.Pz=Qb2.Pz(); file >> B2.Pz;
        B2.Costheta=-costheta.getVal();file >> B2.Costheta;
        B2.Phi=(TMath::Pi()-phi.getVal());file >> B2.Phi;

        K.E=Qk.E();   file >> K.E;
        K.P=Qk.P();   file >> K.P;
        K.Px=Qk.Px(); file >> K.Px;
        K.Py=Qk.Py(); file >> K.Py;
        K.Pz=Qk.Pz(); file >> K.Pz;
        K.Costheta=costheta1.getVal();file >> K.Costheta;
        K.Phi=phi1.getVal();file >> K.Phi;

        Nu1.E=Qnu1.E();   file >> Nu1.E;
        Nu1.P=Qnu1.P();   file >> Nu1.P;
        Nu1.Px=Qnu1.Px(); file >> Nu1.Px;
        Nu1.Py=Qnu1.Py(); file >> Nu1.Py;
        Nu1.Pz=Qnu1.Pz(); file >> Nu1.Pz;
        Nu1.Costheta=costheta2.getVal();file >> Nu1.Costheta;
        Nu1.Phi=phi2.getVal();file >> Nu1.Phi;

        Nu2.E=Qnu2.E();   file >> Nu2.E;
        Nu2.P=Qnu2.P();   file >> Nu2.P;
        Nu2.Px=Qnu2.Px(); file >> Nu2.Px;
        Nu2.Py=Qnu2.Py(); file >> Nu2.Py;
        Nu2.Pz=Qnu2.Pz(); file >> Nu2.Pz;
        Nu2.Costheta=-costheta2.getVal();file >> Nu2.Costheta;
        Nu2.Phi=(TMath::Pi()-phi2.getVal());file >> Nu2.Phi;

        NuNu.E=Qnunu.E();   file >> NuNu.E;
        NuNu.P=Qnunu.P();   file >> NuNu.P;
        NuNu.Px=Qnunu.Px(); file >> NuNu.Px;
        NuNu.Py=Qnunu.Py(); file >> NuNu.Py;
        NuNu.Pz=Qnunu.Pz(); file >> NuNu.Pz;
        NuNu.Costheta=-costheta1.getVal();file >> NuNu.Costheta;
        NuNu.Phi=(TMath::Pi()-phi1.getVal());file >> NuNu.Phi;

        tree->Fill();
        
    }
    file.close();
    tree->Write();
    ftree->Close();
}