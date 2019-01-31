#include <iostream>
#include <TRandom3.h>
#include "fct.cpp"
using namespace std;

RooAbsPdf * m12distrib(RooRealVar& m12)
{
    RooGenericPdf * m12fct = new RooGenericPdf("m12_fct","m12*sqrt((pow(5.27931,2)-pow(m12+0.89166,2))*(pow(5.27931,2)-pow(m12-0.89166,2)))",m12);
    return m12fct;
}

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
    Dalitz dalitz;

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
    tree->Branch("Dalitz", &dalitz, "m12_2/D:m23_2/D");

    //All the constant in our simulation
    RooRealVar Ep("Ep","Positron energy (GeV)",4.0);
    RooRealVar Ee("Ee","Electron energy (GeV)",7.0);
    RooRealVar angle("angle","Angle between e+e- beam (rad)",0.083);
    RooRealVar mb("mb","mass b+/b-",5.27931);
    RooRealVar mk("mk","mass k*+/-",0.89166);
    //RooRealVar mb("mb","mass b0",5.27962);
    //RooRealVar mk("mk","mass k*0",0.89581);

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

    RooRealVar m12("m12","mass 2nu system",0,0,4.38765);

    //Calculate the boost for the B mesons
    TVector3 betaU = Boost(Qu);

    //m12 distribution
    RooAbsPdf* m12distribution = m12distrib(m12);
    RooAbsReal * m12Cdf = m12distribution->createCdf(m12);
    TF1 * m12CdfTH1 = m12Cdf->asTF(m12);

    //List for different resolution comparison
    double_t QnunuE5[nbevts];
    double_t QnunuE10[nbevts];
    double_t QnunuE30[nbevts];
    int k=0;
    double_t sigEmiss;

    for(int i=0;i<nbevts;i++)
    {
        //Get values for all the random variable
        costheta.setVal(random->Uniform(costheta.getMin(),costheta.getMax()));
        phi.setVal(random->Uniform(phi.getMin(),phi.getMax()));

        costheta1.setVal(random->Uniform(costheta1.getMin(),costheta1.getMax()));
        phi1.setVal(random->Uniform(phi1.getMin(),phi1.getMax()));

        costheta2.setVal(random->Uniform(costheta2.getMin(),costheta2.getMax()));
        phi2.setVal(random->Uniform(phi2.getMin(),phi2.getMax()));
        
        double_t u = random->Uniform(0,1);
        
        //Calculate the 4-vector for B & B_
        TLorentzVector Qb1 = TransfoLorentz(betaU,costheta,phi,pb,q/2);
        TLorentzVector Qb2 = TransfoLorentz(betaU,costheta,phi,-pb,q/2);

        //Calculate the 4-vector for K & 2nu
        TVector3 betaB1 = Boost(Qb1);

        m12.setVal(m12CdfTH1->GetX(u));

        double_t pK = sqrt((pow(mb.getVal(),2)-pow(m12.getVal()+mk.getVal(),2))*(pow(mb.getVal(),2)-pow(m12.getVal()-mk.getVal(),2)))/(2*mb.getVal());
        double_t EK = sqrt(mk.getVal()*mk.getVal()+pK*pK);
        
        double_t pnunu = -pK;
        double_t Enunu = mb.getVal()-EK;

        TLorentzVector Qk = TransfoLorentz(betaB1, costheta1, phi1, pK, EK);
        TLorentzVector Qnunu = TransfoLorentz(betaB1, costheta1, phi1, pnunu, Enunu);

        //Resolution study
        sigEmiss = random->Gaus(0,Qnunu.E()*0.05);
        QnunuE5[k]=Qnunu.E()+sigEmiss;
        sigEmiss = random->Gaus(0,Qnunu.E()*0.1);
        QnunuE10[k]=Qnunu.E()+sigEmiss;
        sigEmiss = random->Gaus(0,Qnunu.E()*0.3);
        QnunuE30[k]=Qnunu.E()+sigEmiss;
        k+=1;

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

        dalitz.m12=pow(m12.getVal(),2);file >> dalitz.m12;
        dalitz.m23=(Qk+Qnu2)*(Qk+Qnu2);file >> dalitz.m23;

        tree->Fill();
        
    }

    //Purity/Efficiency/Rejection study
    int Nbin=100;
    double_t xmin=0.0,xmax=10;

    //Different resolution TH1
    TH1D *missE5 = new TH1D("missE5","missE5",Nbin,xmin,xmax);
    TH1D *missE10 = new TH1D("missE10","missE10",Nbin,xmin,xmax);
    TH1D *missE30 = new TH1D("missE30","missE30",Nbin,xmin,xmax);

    for(int i=0;i<nbevts;i+=1)
    {
        missE5->Fill(QnunuE5[i]);
        missE10->Fill(QnunuE10[i]);
        missE30->Fill(QnunuE30[i]);
    }

    //TrueBkg_com to lab frame
    TFile* f1 = new TFile("missingEnergy.root");
    TTree* t1 = (TTree*)f1->Get("tree");

    TH1D *Bkg = new TH1D("Bkg","Bkg",Nbin,xmin,xmax);
    t1->Draw("E>>Bkg","","goff");

    //Draw efficiency histogram
    TH1D *missE = new TH1D("missE","Missing energy",Nbin,xmin,xmax);
    tree->Draw("NuNuevent.E>>missE","","goff");
    
    TH1D *TH1efficiency = new TH1D("TH1efficiency","Efficiency",Nbin,xmin,xmax);
    TH1D *TH1efficiency5 = new TH1D("TH1efficiency5","Efficiency",Nbin,xmin,xmax);
    TH1D *TH1efficiency10 = new TH1D("TH1efficiency10","Efficiency",Nbin,xmin,xmax);
    TH1D *TH1efficiency30 = new TH1D("TH1efficiency30","Efficiency",Nbin,xmin,xmax);

    double_t Vefficiency[Nbin+1];
    double_t Vefficiency5[Nbin+1];
    double_t Vefficiency10[Nbin+1];
    double_t Vefficiency30[Nbin+1];

    Efficiency(Vefficiency,missE,Nbin,xmin,xmax);
    Efficiency(Vefficiency5,missE5,Nbin,xmin,xmax);
    Efficiency(Vefficiency10,missE10,Nbin,xmin,xmax);
    Efficiency(Vefficiency30,missE30,Nbin,xmin,xmax);

    for(int i=0;i<Nbin+1;i+=1)
    {
        TH1efficiency->SetBinContent(i,Vefficiency[i]);
        TH1efficiency5->SetBinContent(i,Vefficiency5[i]);
        TH1efficiency10->SetBinContent(i,Vefficiency10[i]);
        TH1efficiency30->SetBinContent(i,Vefficiency30[i]);
    }
    

    //Draw rejection histogram
    TH1D *TH1rejection = new TH1D("TH1rejection","rejection",Nbin,xmin,xmax);

    double_t Vrejection[Nbin+1];
    Rejection(Vrejection,Bkg,Nbin,xmin,xmax);

    for(int i=0;i<Nbin+1;i+=1)
    {
        TH1rejection->SetBinContent(i,Vrejection[i]);
    }


    //Add bkg with signal
    TH1D *sumbkgsig = new TH1D("sumbkgsig","sumbkgsig",Nbin,xmin,xmax);
    TH1D *sumbkgsig5 = new TH1D("sumbkgsig5","sumbkgsig",Nbin,xmin,xmax);
    TH1D *sumbkgsig10 = new TH1D("sumbkgsig10","sumbkgsig",Nbin,xmin,xmax);
    TH1D *sumbkgsig30 = new TH1D("sumbkgsig30","sumbkgsig",Nbin,xmin,xmax);
    sumbkgsig->Add(Bkg,missE);
    sumbkgsig5->Add(Bkg,missE5);
    sumbkgsig10->Add(Bkg,missE10);
    sumbkgsig30->Add(Bkg,missE30);
	
    //Draw purity histogram
    TH1D *TH1purity = new TH1D("TH1purity","purity",Nbin,xmin,xmax);
    TH1D *TH1purity5 = new TH1D("TH1purity5","purity",Nbin,xmin,xmax);
    TH1D *TH1purity10 = new TH1D("TH1purity10","purity",Nbin,xmin,xmax);
    TH1D *TH1purity30 = new TH1D("TH1purity30","purity",Nbin,xmin,xmax);
    
    double_t Vpurity[Nbin+1];
    double_t Vpurity5[Nbin+1];
    double_t Vpurity10[Nbin+1];
    double_t Vpurity30[Nbin+1];

    Purity(Vpurity,sumbkgsig,missE,Nbin,xmin,xmax);
    Purity(Vpurity5,sumbkgsig5,missE5,Nbin,xmin,xmax);
    Purity(Vpurity10,sumbkgsig10,missE10,Nbin,xmin,xmax);
    Purity(Vpurity30,sumbkgsig30,missE30,Nbin,xmin,xmax);

    for(int i=0;i<Nbin+1;i+=1)
    {
        TH1purity->SetBinContent(i,Vpurity[i]);
        TH1purity5->SetBinContent(i,Vpurity5[i]);
        TH1purity10->SetBinContent(i,Vpurity10[i]);
        TH1purity30->SetBinContent(i,Vpurity30[i]);
    }


    //Significance study
    TH1D *TH1significance = new TH1D("TH1significance","significance",Nbin,xmin,xmax);
    double_t Vsignificance[Nbin];

    double_t cut_bkg,cut_signal;
    double_t fraction=10;

    for(int i=0;i<Nbin+1;i+=1)
    {
        cut_bkg = Bkg->Integral(i,Nbin+1);
        cut_signal = missE->Integral(i,Nbin+1);

        if((sqrt(cut_signal+fraction*cut_bkg))<0.001)
        {
            Vsignificance[i]=0.0;
        }
        else
        {
            Vsignificance[i] = (sqrt(nbevts)*cut_signal)/(sqrt(cut_signal+fraction*cut_bkg));
        }

        TH1significance->SetBinContent(i,Vsignificance[i]);
    }

    //Draw 2D graph
    TGraph *GraphPurEff = new TGraph (Nbin,Vefficiency,Vpurity);
    GraphPurEff->SetTitle("Purity / Efficiency");
    GraphPurEff->GetHistogram()->GetXaxis()->SetTitle("Efficiency");
    GraphPurEff->GetHistogram()->GetYaxis()->SetTitle("Purity"); 
    GraphPurEff->SetMarkerStyle(2);
    TGraph *GraphRejEff = new TGraph (Nbin,Vefficiency,Vrejection);
    GraphRejEff->SetTitle("Rejection / Efficiency");
    GraphRejEff->GetHistogram()->GetXaxis()->SetTitle("Efficiency");
    GraphRejEff->GetHistogram()->GetYaxis()->SetTitle("Rejection");
    GraphRejEff->SetMarkerStyle(2);
    TGraph *GraphPurEff5 = new TGraph (Nbin,Vefficiency5,Vpurity5);
    GraphPurEff5->SetLineColor(2);
    GraphPurEff5->SetMarkerStyle(2);
    GraphPurEff5->SetMarkerColor(2);
    TGraph *GraphRejEff5 = new TGraph (Nbin,Vefficiency5,Vrejection);
    GraphRejEff5->SetLineColor(2);
    GraphRejEff5->SetMarkerStyle(2);
    GraphRejEff5->SetMarkerColor(2);
    TGraph *GraphPurEff10 = new TGraph (Nbin,Vefficiency10,Vpurity10);
    GraphPurEff10->SetLineColor(4);
    GraphPurEff10->SetMarkerStyle(2);
    GraphPurEff10->SetMarkerColor(4);
    TGraph *GraphRejEff10 = new TGraph (Nbin,Vefficiency10,Vrejection);
    GraphRejEff10->SetLineColor(4);
    GraphRejEff10->SetMarkerStyle(2);
    GraphRejEff10->SetMarkerColor(4);
    TGraph *GraphPurEff30 = new TGraph (Nbin,Vefficiency30,Vpurity30);
    GraphPurEff30->SetLineColor(8);
    GraphPurEff30->SetMarkerStyle(2);
    GraphPurEff30->SetMarkerColor(8);
    TGraph *GraphRejEff30 = new TGraph (Nbin,Vefficiency30,Vrejection);
    GraphRejEff30->SetLineColor(8);
    GraphRejEff30->SetMarkerStyle(2);
    GraphRejEff30->SetMarkerColor(8);

    //Create the canvas to plot the results
    missE5->SetLineColor(2);
    missE10->SetLineColor(4);
    missE30->SetLineColor(8);

    TH1efficiency5->SetLineColor(2);
    TH1efficiency10->SetLineColor(4);
    TH1efficiency30->SetLineColor(8);

    TH1purity->SetMaximum(1.);
    TH1purity5->SetLineColor(2);
    TH1purity10->SetLineColor(4);
    TH1purity30->SetLineColor(8);

    TCanvas* c = new TCanvas("ResultsTIPP","ResultsTIPP",900,450) ;
	c->Divide(2,1);
	c->cd(1);
	gPad->SetLeftMargin(0.15) ;  
    GraphPurEff30->Draw("APL") ;
    GraphPurEff->Draw("same,PL") ;
    GraphPurEff5->Draw("same,PL") ;
    GraphPurEff10->Draw("same,PL") ;
    c->cd(2);
	gPad->SetLeftMargin(0.15) ;  
    GraphRejEff->Draw("APL") ; 
    GraphRejEff5->Draw("same,PL") ;
    GraphRejEff10->Draw("same,PL") ;
    GraphRejEff30->Draw("same,PL") ;

    TCanvas* d = new TCanvas("ResultsTIPP2","ResultsTIPP2",1350,900) ;
	d->Divide(3,2);
	d->cd(1);
	gPad->SetLeftMargin(0.15) ;  
    TH1efficiency->Draw() ;
    TH1efficiency5->Draw("same") ;
    TH1efficiency10->Draw("same") ;
    TH1efficiency30->Draw("same") ;

    d->cd(2);
	gPad->SetLeftMargin(0.15) ;  
    TH1rejection->Draw() ;

    d->cd(3);
	gPad->SetLeftMargin(0.15) ;
    TH1purity->Draw() ;  
    TH1purity30->Draw("same") ;  
    TH1purity5->Draw("same") ;  
    TH1purity10->Draw("same") ;

    d->cd(4);
	gPad->SetLeftMargin(0.15) ;  
    Bkg->Draw() ;

    d->cd(5);
	gPad->SetLeftMargin(0.15) ;  
    missE->Draw() ; 
    missE5->Draw("same"); 
    missE10->Draw("same"); 
    missE30->Draw("same");

    d->cd(6);
	//gPad->SetLeftMargin(0.15) ;  
   // sumbkgsig30->Draw() ;
   // missE30->Draw("same,HIST");
   // gPad->SetLogy();
    TH1significance->Draw();

    TLegend *leg = new TLegend( .58, .65, .9, .9, "Emiss Resolution");
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->AddEntry( GraphPurEff, "0\%", "lp");
    leg->AddEntry( GraphPurEff5, "5\%", "lp");
    leg->AddEntry( GraphPurEff10, "10\%", "lp");
    leg->AddEntry( GraphPurEff30, "30\%", "lp");
    c->cd(1);leg->Draw();
    d->cd(2);leg->Draw();

    //Useful for different signal/bkg comparison
    //ofstream fileSignal("SignalBkgProp.txt");
    //ofstream fileSignal("SignalBkgProp.txt",std::ios_base::app); //without erasing
    //for (int i=0;i<Nbin+1;i++)
    //{
    //    fileSignal<<Vefficiency[i]<<" "<<Vpurity[i]<<" "<<Vrejection[i]<<endl;
    //}
    //fileSignal.close();

    GraphEffSignal(Nbin+1);

    tree->Fill();
    ftree->Write();

    //Draw the Dalitz plot
    //tree->Draw("m23_2:m12_2>>(250,0,25,300,0,30)","","colz");

    //ftree->Close();
    //f1->Close();
    //file.Close();
}