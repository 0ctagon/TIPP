#include <iostream>
#include <TRandom3.h>
#include "fct.cpp"
using namespace std;

//m12 distribution (function g(mnunu) in the article)
RooAbsPdf * m12distrib(RooRealVar& m12)
{
    RooGenericPdf * m12fct = new RooGenericPdf("m12_fct","m12*sqrt((pow(5.27962,2)-pow(m12+0.89581,2))*(pow(5.27962,2)-pow(m12-0.89581,2)))",m12);
    return m12fct;
}


void MC()
{
    //Set seed for the random variables
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
    Dalitz dalitz;

    //Number of event in the background and in our simulation
    int nbkg = 722287;
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
    //This branch is used to draw the Dalitz plot of the B->Knunu decay
    tree->Branch("Dalitz", &dalitz, "m12_2/D:m23_2/D");


    //All the constant in our simulation
    RooRealVar Ep("Ep","Positron energy (GeV)",4.0);
    RooRealVar Ee("Ee","Electron energy (GeV)",7.0);
    RooRealVar angle("angle","Angle between e+e- beam (rad)",0.083);
    //RooRealVar mb("mb","mass b+/b-(GeV)",5.27931);
    //RooRealVar mk("mk","mass k*+/-(GeV)",0.89166);
    RooRealVar mb("mb","mass B0 (GeV)",5.27962);
    RooRealVar mk("mk","mass k*0 (GeV)",0.89581);


    //Define the 4-vectors Qe/Qp for e-/e+ & Qu for Upsilon 4S
    TLorentzVector Qp(-Ep.getVal(),0,0,Ep.getVal());
    TLorentzVector Qe(Ee.getVal()*cos(angle.getVal()),-Ee.getVal()*sin(angle.getVal()),0,Ee.getVal());
    TLorentzVector Qu = Qp + Qe;


    //Calculate q = sqrt(s) and pb the impulsion of the B mesons
    double_t q = sqrt(Qu*Qu);
    double_t pb = sqrt(q*q-4*mb.getVal()*mb.getVal())/2;


    //Definition of all the random variables used in the simulation
    RooRealVar costheta("costheta","angle between B-/B in the COM",0,-1,1);
    RooRealVar phi("phi","angle between B-/B in the COM",0,0,TMath::TwoPi());

    RooRealVar costheta1("costheta1","angle between K/2nu in the COM of B",0,-1,1);
    RooRealVar phi1("phi","angle between K/2nu in the COM of B",0,0,TMath::TwoPi());

    RooRealVar costheta2("costheta2","angle between 2nu in the COM of B",0,-1,1);
    RooRealVar phi2("phi2","angle between 2nu in the COM of B",0,0,TMath::TwoPi());

    RooRealVar m12("m12","mass 2nu system",0,0,4.38381);//Defined between 0 and mb-mk=4.38381

    //Calculate the boost for the B mesons
    TVector3 betaU = Boost(Qu);

    //m12 distribution (g in the article)
    RooAbsPdf* m12distribution = m12distrib(m12);
    //We need the Cumulative distribution function m12cdf because we use the Inverse transform method to randomly distribute m12 following m12distrib
    RooAbsReal * m12Cdf = m12distribution->createCdf(m12);
    TF1 * m12CdfTH1 = m12Cdf->asTF(m12);

    //Define list for different resolution on Emiss
    double_t* QnunuE10 = new double_t[nbevts];
    double_t* QnunuE20 = new double_t[nbevts];
    double_t* QnunuE30 = new double_t[nbevts];
    int k=0;
    double_t sigEmiss;


    for(int i=0;i<nbevts;i++)
    {
        ////Get values for all the random variable
        costheta.setVal(random->Uniform(costheta.getMin(),costheta.getMax()));
        phi.setVal(random->Uniform(phi.getMin(),phi.getMax()));

        costheta1.setVal(random->Uniform(costheta1.getMin(),costheta1.getMax()));
        phi1.setVal(random->Uniform(phi1.getMin(),phi1.getMax()));

        costheta2.setVal(random->Uniform(costheta2.getMin(),costheta2.getMax()));
        phi2.setVal(random->Uniform(phi2.getMin(),phi2.getMax()));
        
        //For m12 we use the Inverse transform method
        double_t u = random->Uniform(0,1);
        

        ////Calculate the 4-vector for B & Bbar
        TLorentzVector Qb1 = TransfoLorentz(betaU,costheta,phi,pb,q/2);
        TLorentzVector Qb2 = TransfoLorentz(betaU,costheta,phi,-pb,q/2);


        ////Calculate the 4-vector for K (Qk) & 2nu system (Qnunu)
        //Boost of the B
        TVector3 betaB1 = Boost(Qb1);

        //Get the value for m12
        m12.setVal(m12CdfTH1->GetX(u));

        //Calculate the impulsion and energy of K
        double_t pK = sqrt((pow(mb.getVal(),2)-pow(m12.getVal()+mk.getVal(),2))*(pow(mb.getVal(),2)-pow(m12.getVal()-mk.getVal(),2)))/(2*mb.getVal());
        double_t EK = sqrt(mk.getVal()*mk.getVal()+pK*pK);
        
        ////Calculate the impulsion and energy of 2nu system
        double_t pnunu = -pK;
        double_t Enunu = mb.getVal()-EK;

        TLorentzVector Qk = TransfoLorentz(betaB1, costheta1, phi1, pK, EK);
        TLorentzVector Qnunu = TransfoLorentz(betaB1, costheta1, phi1, pnunu, Enunu);


        ////Resolution study on Emiss
        //10% resolution
        sigEmiss = random->Gaus(0,Qnunu.E()*0.1);
        QnunuE10[k]=Qnunu.E()+sigEmiss;
        //20% resolution
        sigEmiss = random->Gaus(0,Qnunu.E()*0.2);
        QnunuE20[k]=Qnunu.E()+sigEmiss;
        //30% resolution
        sigEmiss = random->Gaus(0,Qnunu.E()*0.3);
        QnunuE30[k]=Qnunu.E()+sigEmiss;
        k+=1;


        //Calculate the 4-vector for nu1 & nu2 (not used)
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

        dalitz.m12=pow(m12.getVal(),2);file >> dalitz.m12;
        dalitz.m23=(Qk+Qnu2)*(Qk+Qnu2);file >> dalitz.m23;

        tree->Fill();
        
    }

    ////Purity/Efficiency/Rejection study

    int Nbin=100;
    double_t xmin=0.0,xmax=7.;

    //TH1 for differents resolutions
    TH1D *missE10 = new TH1D("missE10","missE10",Nbin,xmin,xmax);
    TH1D *missE20 = new TH1D("missE20","missE20",Nbin,xmin,xmax);
    TH1D *missE30 = new TH1D("missE30","missE30",Nbin,xmin,xmax);

    for(int i=0;i<nbevts;i+=1)
    {
        missE10->Fill(QnunuE10[i]);
        missE20->Fill(QnunuE20[i]);
        missE30->Fill(QnunuE30[i]);
    }

    //Get the total Energy background Bkg from the .root file
    TFile* f1 = new TFile("missingEnergy.root");
    TTree* t1 = (TTree*)f1->Get("tree");

    TH1D *Bkg = new TH1D("Bkg","Bkg",Nbin,xmin,xmax);
    t1->Draw("E>>Bkg","","goff");


    //Draw efficiency histogram for differents resolutions
    TH1D *missE = new TH1D("missE","Missing energy",Nbin,xmin,xmax);
    tree->Draw("NuNuevent.E>>missE","","goff");
    
    TH1D *TH1efficiency = new TH1D("TH1efficiency","Efficiency",Nbin,xmin,xmax);
    TH1D *TH1efficiency10 = new TH1D("TH1efficiency10","Efficiency",Nbin,xmin,xmax);
    TH1D *TH1efficiency20 = new TH1D("TH1efficiency20","Efficiency",Nbin,xmin,xmax);
    TH1D *TH1efficiency30 = new TH1D("TH1efficiency30","Efficiency",Nbin,xmin,xmax);

    double_t Vefficiency[Nbin+1];
    double_t Vefficiency10[Nbin+1];
    double_t Vefficiency20[Nbin+1];
    double_t Vefficiency30[Nbin+1];

    Efficiency(Vefficiency,TH1efficiency,missE,Nbin,xmin,xmax);
    Efficiency(Vefficiency10,TH1efficiency10,missE10,Nbin,xmin,xmax);
    Efficiency(Vefficiency20,TH1efficiency20,missE20,Nbin,xmin,xmax);
    Efficiency(Vefficiency30,TH1efficiency30,missE30,Nbin,xmin,xmax);
    

    //Draw rejection histogram (doesn't change for differents resolutions)
    TH1D *TH1rejection = new TH1D("TH1rejection","rejection",Nbin,xmin,xmax);

    double_t Vrejection[Nbin+1];
    Rejection(Vrejection,TH1rejection,Bkg,Nbin,xmin,xmax);


    //Add bkg with signal for differents resolutions
    TH1D *sumbkgsig = new TH1D("sumbkgsig","sumbkgsig",Nbin,xmin,xmax);
    TH1D *sumbkgsig10 = new TH1D("sumbkgsig10","sumbkgsig",Nbin,xmin,xmax);
    TH1D *sumbkgsig20 = new TH1D("sumbkgsig20","sumbkgsig",Nbin,xmin,xmax);
    TH1D *sumbkgsig30 = new TH1D("sumbkgsig30","sumbkgsig",Nbin,xmin,xmax);
    sumbkgsig->Add(Bkg,missE);
    sumbkgsig10->Add(Bkg,missE10);
    sumbkgsig20->Add(Bkg,missE20);
    sumbkgsig30->Add(Bkg,missE30);
	

    //Draw purity histogram for differents resolutions
    TH1D *TH1purity = new TH1D("TH1purity","purity",Nbin,xmin,xmax);
    TH1D *TH1purity10 = new TH1D("TH1purity10","purity",Nbin,xmin,xmax);
    TH1D *TH1purity20 = new TH1D("TH1purity20","purity",Nbin,xmin,xmax);
    TH1D *TH1purity30 = new TH1D("TH1purity30","purity",Nbin,xmin,xmax);
    
    double_t Vpurity[Nbin+1];
    double_t Vpurity10[Nbin+1];
    double_t Vpurity20[Nbin+1];
    double_t Vpurity30[Nbin+1];

    Purity(Vpurity,TH1purity,sumbkgsig,missE,Nbin,xmin,xmax);
    Purity(Vpurity10,TH1purity10,sumbkgsig10,missE10,Nbin,xmin,xmax);
    Purity(Vpurity20,TH1purity20,sumbkgsig20,missE20,Nbin,xmin,xmax);
    Purity(Vpurity30,TH1purity30,sumbkgsig30,missE30,Nbin,xmin,xmax);


    ////Significance study
    TH1D *TH1significance1 = new TH1D("TH1significance1","purity",Nbin,xmin,xmax);
    TH1D *TH1significance2 = new TH1D("TH1significance2","purity",Nbin,xmin,xmax);
    TH1D *TH1significance3 = new TH1D("TH1significance3","purity",Nbin,xmin,xmax);
    TH1D *TH1significance4 = new TH1D("TH1significance4","purity",Nbin,xmin,xmax);
    TH1D *TH1significance5 = new TH1D("TH1significance5","purity",Nbin,xmin,xmax);

    //table with significance
    double_t S10f1 = SignificanceTH(TH1significance1,TH1efficiency10,TH1rejection,Nbin,nbkg,1,xmin,xmax,false);
    double_t S10f10 = SignificanceTH(TH1significance2,TH1efficiency10,TH1rejection,Nbin,nbkg,10,xmin,xmax,false);
    double_t S10f100 = SignificanceTH(TH1significance3,TH1efficiency10,TH1rejection,Nbin,nbkg,100,xmin,xmax,false);
    double_t S10f1000 = SignificanceTH(TH1significance4,TH1efficiency10,TH1rejection,Nbin,nbkg,1000,xmin,xmax,false);
    double_t S10f100000 = SignificanceTH(TH1significance5,TH1efficiency10,TH1rejection,Nbin,nbkg,100000,xmin,xmax,false);
    
    double_t S20f1 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,1,xmin,xmax,false);
    double_t S20f10 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,10,xmin,xmax,false);
    double_t S20f100 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,100,xmin,xmax,false);
    double_t S20f1000 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,1000,xmin,xmax,false);
    double_t S20f100000 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,100000,xmin,xmax,false);

    double_t S30f1 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,1,xmin,xmax,false);
    double_t S30f10 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,10,xmin,xmax,false);
    double_t S30f100 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,100,xmin,xmax,false);
    double_t S30f1000 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,1000,xmin,xmax,false);
    double_t S30f100000 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,100000,xmin,xmax,false);

    cout<<"____________________"<<endl;
    cout<<"Significance and cut values for different background to signal fraction and Resolution on Emiss"<<endl;
    cout<<endl;
    cout<<"Significance values"<<endl;
    cout<<"Fraction:        "<<"1"<<"\t"<<"10"<<"\t"<<"100"<<"\t"<<"1000"<<"\t"<<"100000"<<endl;
    cout<<endl;
    cout<<"Resolution: 10\%  "<<setprecision(3)<<S10f1<<"\t"<<S10f10<<"\t"<<S10f100<<"\t"<<S10f1000<<"\t"<<S10f100000<<endl;
    cout<<"Resolution: 20\%  "<<setprecision(3)<<S20f1<<"\t"<<S20f10<<"\t"<<S20f100<<"\t"<<S20f1000<<"\t"<<S20f100000<<endl;
    cout<<"Resolution: 30\%  "<<setprecision(3)<<S30f1<<"\t"<<S30f10<<"\t"<<S30f100<<"\t"<<S30f1000<<"\t"<<S30f100000<<endl; 
    cout<<endl;

    //table with cut
    S10f1 = Significance(TH1efficiency10,TH1rejection,Nbin,nbkg,1,xmin,xmax,true);
    S10f10 = Significance(TH1efficiency10,TH1rejection,Nbin,nbkg,10,xmin,xmax,true);
    S10f100 = Significance(TH1efficiency10,TH1rejection,Nbin,nbkg,100,xmin,xmax,true);
    S10f1000 = Significance(TH1efficiency10,TH1rejection,Nbin,nbkg,1000,xmin,xmax,true);
    S10f100000 = Significance(TH1efficiency10,TH1rejection,Nbin,nbkg,100000,xmin,xmax,true);
    
    S20f1 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,1,xmin,xmax,true);
    S20f10 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,10,xmin,xmax,true);
    S20f100 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,100,xmin,xmax,true);
    S20f1000 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,1000,xmin,xmax,true);
    S20f100000 = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,100000,xmin,xmax,true);

    S30f1 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,1,xmin,xmax,true);
    S30f10 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,10,xmin,xmax,true);
    S30f100 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,100,xmin,xmax,true);
    S30f1000 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,1000,xmin,xmax,true);
    S30f100000 = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,100000,xmin,xmax,true);

    cout<<endl;
    cout<<"Cut values"<<endl;
    cout<<"Fraction:        "<<"1"<<"\t"<<"10"<<"\t"<<"100"<<"\t"<<"1000"<<"\t"<<"100000"<<endl;
    cout<<endl;
    cout<<"Resolution: 10\%  "<<setprecision(3)<<S10f1<<"\t"<<S10f10<<"\t"<<S10f100<<"\t"<<S10f1000<<"\t"<<S10f100000<<endl;
    cout<<"Resolution: 20\%  "<<setprecision(3)<<S20f1<<"\t"<<S20f10<<"\t"<<S20f100<<"\t"<<S20f1000<<"\t"<<S20f100000<<endl;
    cout<<"Resolution: 30\%  "<<setprecision(3)<<S30f1<<"\t"<<S30f10<<"\t"<<S30f100<<"\t"<<S30f1000<<"\t"<<S30f100000<<endl; 
    cout<<endl;


    ////#bkg to #signal fraction needed to get 3 and 5sigma
    cout<<"____________________"<<endl;
    cout<<"Background to signal fraction and significance study"<<endl;
	cout<<endl;
	cout<<"10\% Resolution"<<endl;
    double_t S;
	for(int i=1000;i>=100;i--)
	{
		S = Significance(TH1efficiency10,TH1rejection,Nbin,nbkg,i,xmin,xmax,false);
        if((S<=3.01) && (S>=2.99))
		{
			cout<<"Fraction = "<<i<<" : "<<"Significance = "<<S<<endl;
		}
		if((S<=5.02) && (S>=4.98))
		{
			cout<<"Fraction = "<<i<<" : "<<"Significance = "<<S<<endl;
		}
	}
	cout<<endl;
	cout<<"20\% Resolution"<<endl;
	for(int i=1000;i>=100;i--)
	{
		S = Significance(TH1efficiency20,TH1rejection,Nbin,nbkg,i,xmin,xmax,false);
		if((S<=3.01) && (S>=2.99))
		{
			cout<<"Fraction = "<<i<<" : "<<"Significance = "<<S<<endl;
		}
		if((S<=5.02) && (S>=4.98))
		{
			cout<<"Fraction = "<<i<<" : "<<"Significance = "<<S<<endl;
		}
	}
	cout<<endl;
	cout<<"30\% Resolution"<<endl;
	for(int i=1000;i>=100;i--)
	{
		S = Significance(TH1efficiency30,TH1rejection,Nbin,nbkg,i,xmin,xmax,false);
		if((S<=3.01) && (S>=2.99))
		{
			cout<<"Fraction = "<<i<<" : "<<"Significance = "<<S<<endl;
		}
		if((S<=5.02) && (S>=4.98))
		{
			cout<<"Fraction = "<<i<<" : "<<"Significance = "<<S<<endl;
		}
	}
    cout<<endl;


    //Draw 2D graph
    TGraph *GraphPurEff = new TGraph (Nbin,Vefficiency,Vpurity);
    GraphPurEff->SetMarkerStyle(2);
    TGraph *GraphRejEff = new TGraph (Nbin,Vefficiency,Vrejection);
    GraphRejEff->SetTitle("Rejection / Efficiency");
    GraphRejEff->GetHistogram()->GetXaxis()->SetTitle("Signal Efficiency");
    GraphRejEff->GetHistogram()->GetYaxis()->SetTitle("Bkg Rejection");
    GraphRejEff->SetMarkerStyle(2);
    TGraph *GraphPurEff10 = new TGraph (Nbin,Vefficiency10,Vpurity10);
    GraphPurEff10->SetLineColor(2);
    GraphPurEff10->SetMarkerStyle(2);
    GraphPurEff10->SetMarkerColor(2);
    TGraph *GraphRejEff10 = new TGraph (Nbin,Vefficiency10,Vrejection);
    GraphRejEff10->SetLineColor(2);
    GraphRejEff10->SetMarkerStyle(2);
    GraphRejEff10->SetMarkerColor(2);
    TGraph *GraphPurEff20 = new TGraph (Nbin,Vefficiency20,Vpurity20);
    GraphPurEff20->SetLineColor(4);
    GraphPurEff20->SetMarkerStyle(2);
    GraphPurEff20->SetMarkerColor(4);
    TGraph *GraphRejEff20 = new TGraph (Nbin,Vefficiency20,Vrejection);
    GraphRejEff20->SetLineColor(4);
    GraphRejEff20->SetMarkerStyle(2);
    GraphRejEff20->SetMarkerColor(4);
    TGraph *GraphPurEff30 = new TGraph (Nbin,Vefficiency30,Vpurity30);
    GraphPurEff30->SetTitle("Purity / Efficiency");
    GraphPurEff30->GetHistogram()->GetXaxis()->SetTitle("Signal Efficiency");
    GraphPurEff30->GetHistogram()->GetYaxis()->SetTitle("Purity"); 
    GraphPurEff30->SetLineColor(8);
    GraphPurEff30->SetMarkerStyle(2);
    GraphPurEff30->SetMarkerColor(8);
    TGraph *GraphRejEff30 = new TGraph (Nbin,Vefficiency30,Vrejection);
    GraphRejEff30->SetLineColor(8);
    GraphRejEff30->SetMarkerStyle(2);
    GraphRejEff30->SetMarkerColor(8);

    //Create the canvas to plot the results
    missE10->SetLineColor(2);
    missE20->SetLineColor(4);
    missE30->SetLineColor(8);

    TH1efficiency->SetTitle("Efficiency");
    TH1efficiency->GetXaxis()->SetTitle("Cut MissE [GeV]");
    TH1efficiency->GetYaxis()->SetTitle("Efficiency");
    TH1efficiency10->SetLineColor(2);
    TH1efficiency20->SetLineColor(4);
    TH1efficiency30->SetLineColor(8);

    TH1purity->SetTitle("Purity");
    TH1purity->GetXaxis()->SetTitle("Cut MissE [GeV]");
    TH1purity->GetYaxis()->SetTitle("Purity");
    TH1purity->SetMaximum(1.);
    TH1purity10->SetLineColor(2);
    TH1purity20->SetLineColor(4);
    TH1purity30->SetLineColor(8);

    TH1rejection->SetTitle("Rejection");
    TH1rejection->GetXaxis()->SetTitle("Cut Bkg [GeV]");
    TH1rejection->GetYaxis()->SetTitle("Rejection");

    Bkg->SetTitle("Background");
    Bkg->GetXaxis()->SetTitle("E [GeV]");
    Bkg->GetYaxis()->SetTitle("Events/0.07 GeV");

    missE->SetTitle("Missing Energy");
    missE->GetXaxis()->SetTitle("Enunu [GeV]");
    missE->GetYaxis()->SetTitle("Events/0.07 GeV");

    TH1significance1->SetTitle("Significance (10\% Resolution)");
    TH1significance1->GetXaxis()->SetTitle("Cut [GeV]");
    TH1significance1->GetYaxis()->SetTitle("Significance");
    TH1significance1->SetMarkerStyle(0);
    TH1significance2->SetLineColor(2);
    TH1significance2->SetMarkerStyle(0);
    TH1significance3->SetLineColor(4);
    TH1significance3->SetMarkerStyle(0);
    TH1significance4->SetLineColor(8);
    TH1significance4->SetMarkerStyle(0);
    TH1significance5->SetLineColor(28);
    TH1significance5->SetMarkerStyle(0);

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("ResultsTIPP","ResultsTIPP",900,450) ;
	c->Divide(2,1);
	c->cd(1);
	gPad->SetLeftMargin(0.15) ;  
    GraphPurEff30->Draw("APL") ;
    GraphPurEff->Draw("same,PL") ;
    GraphPurEff10->Draw("same,PL") ;
    GraphPurEff20->Draw("same,PL") ;
    c->cd(2);
	gPad->SetLeftMargin(0.15) ;  
    GraphRejEff->Draw("APL") ; 
    GraphRejEff10->Draw("same,PL") ;
    GraphRejEff20->Draw("same,PL") ;
    GraphRejEff30->Draw("same,PL") ;

    TCanvas* d = new TCanvas("ResultsTIPP2","ResultsTIPP2",1350,900) ;
	d->Divide(3,2);
	d->cd(1);
	gPad->SetLeftMargin(0.15) ;  
    TH1efficiency->Draw() ;
    TH1efficiency10->Draw("same") ;
    TH1efficiency20->Draw("same") ;
    TH1efficiency30->Draw("same") ;

    d->cd(2);
	gPad->SetLeftMargin(0.15) ;  
    TH1rejection->Draw() ;

    d->cd(3);
	gPad->SetLeftMargin(0.15) ;
    TH1purity->Draw() ;  
    TH1purity30->Draw("same") ;  
    TH1purity10->Draw("same") ;  
    TH1purity20->Draw("same") ;

    d->cd(4);
	gPad->SetLeftMargin(0.15) ;  
    Bkg->Draw() ;

    d->cd(5);
	gPad->SetLeftMargin(0.15) ;  
    missE->Draw() ; 
    missE10->Draw("same"); 
    missE20->Draw("same"); 
    missE30->Draw("same");

    d->cd(6);
	gPad->SetLeftMargin(0.15) ; 
    gPad->SetLogy();
    TH1significance1->Draw() ;
    TH1significance2->Draw("same") ;
    TH1significance3->Draw("same") ;
    TH1significance4->Draw("same") ;
    TH1significance5->Draw("same") ;

    //Add legend
    TLegend *leg = new TLegend( .58, .65, .9, .9, "Emiss Resolution");
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->AddEntry( GraphPurEff, "0\%", "lp");
    leg->AddEntry( GraphPurEff10, "10\%", "lp");
    leg->AddEntry( GraphPurEff20, "20\%", "lp");
    leg->AddEntry( GraphPurEff30, "30\%", "lp");
    c->cd(1);leg->Draw();
    c->cd(2);leg->Draw();
    d->cd(1);leg->Draw();
    d->cd(3);leg->Draw();
    d->cd(5);leg->Draw();

    TLegend *leg1 = new TLegend( .58, .65, .9, .9, "#Bkg/#Signal");
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1001);
    leg1->AddEntry( TH1significance1, "1", "lp");
    leg1->AddEntry( TH1significance2, "10", "lp");
    leg1->AddEntry( TH1significance3, "100", "lp");
    leg1->AddEntry( TH1significance4, "1000", "lp");
    leg1->AddEntry( TH1significance5, "100000", "lp");
    d->cd(6);leg1->Draw();

    tree->Fill();
    ftree->Write();
}