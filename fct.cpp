#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include <iostream>
using namespace std;

TVector3 Boost(TLorentzVector A)
{
    TVector3 beta(A.Px()/A.E(),A.Py()/A.E(),A.Pz()/A.E());
    return beta;
}

TLorentzVector TransfoLorentz(TVector3 b, RooRealVar costheta, RooRealVar phi, double_t p, double_t q)
{
    double_t _Px = p*costheta.getVal();
    double_t _Py = p*(sqrt(1-costheta.getVal()*costheta.getVal()))*cos(phi.getVal());
    double_t _Pz = p*(sqrt(1-costheta.getVal()*costheta.getVal()))*sin(phi.getVal());
    double_t _E = q;

    double_t gamma = 1/sqrt(1-b*b);
    double_t beta2 = b*b;

    double_t p_l = gamma*b.X()*_E+(gamma-1)*((b.X()*b.Y())/beta2)*_Py+(gamma-1)*((b.X()*b.Z())/beta2)*_Pz+(1+(gamma-1)*((b.X()*b.X())/beta2))*_Px;
    double_t p_t = gamma*b.Y()*_E+(gamma-1)*((b.Y()*b.X())/beta2)*_Px+(gamma-1)*((b.Z()*b.Y())/beta2)*_Pz+(1+(gamma-1)*((b.Y()*b.Y())/beta2))*_Py;
    double_t p_z = gamma*b.Z()*_E+(gamma-1)*((b.Z()*b.X())/beta2)*_Px+(gamma-1)*((b.Z()*b.Y())/beta2)*_Py+(1+(gamma-1)*((b.Z()*b.Z())/beta2))*_Pz;
    double_t E = gamma*(_E+b.X()*_Px+b.Y()*_Py+b.Z()*_Pz);

    TLorentzVector Q(p_l,p_t,p_z,E);
    return Q;
}


void GraphEffSignal(int N)
{
    double_t Veff1[N],Vpur1[N],Vrej1[N];
    double_t Veff2[N],Vpur2[N],Vrej2[N];
    double_t Veff3[N],Vpur3[N],Vrej3[N];
    double_t Veff4[N],Vpur4[N],Vrej4[N];
    double_t Veff5[N],Vpur5[N],Vrej5[N];

    ifstream Vect("SignalBkgProp.txt");
    for (int i=0;i<N;i++)
    {
        Vect >> Veff1[i] >> Vpur1[i] >> Vrej1[i];
    }
    for (int i=0;i<N;i++)
    {
        Vect >> Veff2[i] >> Vpur2[i] >> Vrej2[i];
    }
    for (int i=0;i<N;i++)
    {
        Vect >> Veff3[i] >> Vpur3[i] >> Vrej3[i];
    }
    for (int i=0;i<N;i++)
    {
        Vect >> Veff4[i] >> Vpur4[i] >> Vrej4[i];
    }
    for (int i=0;i<N;i++)
    {
        Vect >> Veff5[i] >> Vpur5[i] >> Vrej5[i];
    }

    TGraph *GraphPurEff1 = new TGraph (N,Veff1,Vpur1);
    GraphPurEff1->SetTitle("Purity / Efficiency");
    GraphPurEff1->GetHistogram()->GetXaxis()->SetTitle("Efficiency");
    GraphPurEff1->GetHistogram()->GetYaxis()->SetTitle("Purity");
    GraphPurEff1->GetHistogram()->SetMaximum(1.);
    GraphPurEff1->SetLineColor(2);
    GraphPurEff1->SetMarkerStyle(2);
    GraphPurEff1->SetMarkerColor(2);
    TGraph *GraphRejEff1 = new TGraph (N,Veff1,Vrej1);
    GraphRejEff1->SetTitle("Rejection / Efficiency");
    GraphRejEff1->GetHistogram()->GetXaxis()->SetTitle("Efficiency");
    GraphRejEff1->GetHistogram()->GetYaxis()->SetTitle("Rejection");
    GraphRejEff1->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
    GraphRejEff1->SetLineColor(2);
    GraphRejEff1->SetMarkerStyle(2);
    GraphRejEff1->SetMarkerColor(2);

    TGraph *GraphPurEff2 = new TGraph (N,Veff2,Vpur2);
    GraphPurEff2->SetLineColor(4);
    GraphPurEff2->SetMarkerStyle(2);
    GraphPurEff2->SetMarkerColor(4);
    TGraph *GraphRejEff2 = new TGraph (N,Veff2,Vrej2);
    GraphRejEff2->SetLineColor(4);
    GraphRejEff2->SetMarkerStyle(2);
    GraphRejEff2->SetMarkerColor(4);
    TGraph *GraphPurEff3 = new TGraph (N,Veff3,Vpur3);
    GraphPurEff3->SetLineColor(8);
    GraphPurEff3->SetMarkerStyle(2);
    GraphPurEff3->SetMarkerColor(8);
    TGraph *GraphRejEff3 = new TGraph (N,Veff3,Vrej3);
    GraphRejEff3->SetLineColor(8);
    GraphRejEff3->SetMarkerStyle(2);
    GraphRejEff3->SetMarkerColor(8);
    TGraph *GraphPurEff4 = new TGraph (N,Veff4,Vpur4);
    GraphPurEff4->SetLineColor(1);
    GraphPurEff4->SetMarkerStyle(2);
    GraphPurEff4->SetMarkerColor(1);
    TGraph *GraphRejEff4 = new TGraph (N,Veff4,Vrej4);
    GraphRejEff4->SetLineColor(1);
    GraphRejEff4->SetMarkerStyle(2);
    GraphRejEff4->SetMarkerColor(1);
    TGraph *GraphPurEff5 = new TGraph (N,Veff5,Vpur5);
    GraphPurEff5->SetLineColor(28);
    GraphPurEff5->SetMarkerStyle(2);
    GraphPurEff5->SetMarkerColor(28);
    TGraph *GraphRejEff5 = new TGraph (N,Veff5,Vrej5);
    GraphRejEff5->SetLineColor(28);
    GraphRejEff5->SetMarkerStyle(2);
    GraphRejEff5->SetMarkerColor(28);

    TCanvas* c = new TCanvas("comparisonSignalBkg","comparison SignalBkg",900,450) ;
	c->Divide(2,1);
	c->cd(1);
	GraphPurEff1->Draw("APL") ;
    GraphPurEff2->Draw("same,PL") ;
    GraphPurEff3->Draw("same,PL") ;
    GraphPurEff4->Draw("same,PL") ;
    //GraphPurEff5->Draw("same,PL") ;
    gPad->SetLogy();
    c->cd(2);
	GraphRejEff1->Draw("APL") ;
    GraphRejEff2->Draw("same,PL") ;
    GraphRejEff3->Draw("same,PL") ;
    GraphRejEff4->Draw("same,PL") ;
    //GraphRejEff5->Draw("same,PL") ;
    TLegend *leg = new TLegend( .58, .65, .9, .9, "#Signal/#Bkg");
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->AddEntry( GraphPurEff1, "1", "lp");
    leg->AddEntry( GraphPurEff2, "10", "lp");
    leg->AddEntry( GraphPurEff3, "100", "lp");
    leg->AddEntry( GraphPurEff4, "1000", "lp");
    //leg->AddEntry( GraphPurEff5, "100000", "lp");
    leg->Draw();
    Vect.close();
}

void Efficiency(double_t* Veff, TH1* hist,int Nbin, double_t xmin, double_t xmax)
{
    double_t full_integral = hist->Integral(), cut_integral;
    double h = xmax/Nbin;

    for(int i=0;i<Nbin+1;i+=1)
    {
        cut_integral = hist->Integral(i,Nbin)/full_integral;
        Veff[i]=cut_integral;
    }
}

void Rejection(double_t* Vrej, TH1* hist, int Nbin, double_t xmin, double_t xmax)
{
    double_t cut_integral;
    double_t full_integral = hist->Integral();

    for(int i=0;i<Nbin+1;i+=1)
    {
        cut_integral = hist->Integral(i,Nbin)/full_integral;
        Vrej[i]=1-cut_integral;
    }
}

void Purity(double_t* Vpur, TH1* histbkg, TH1* histE, int Nbin, double_t xmin, double_t xmax)
{
    double_t sum_integral, cut_integral;

    for(int i=0;i<Nbin+1;i+=1)
    {
        sum_integral = histbkg->Integral(i,Nbin+1);
        cut_integral = histE->Integral(i,Nbin+1);

        if(sum_integral<0.001)
        {
            Vpur[i]=0.0;
        }
        else
        {
            Vpur[i]=cut_integral/sum_integral;
        }
    }
}

struct Particle{
   	double_t E;
    double_t P;
    double_t Px;
    double_t Py;
    double_t Pz;
    double_t Costheta;
    double_t Phi;
};

struct Dalitz{
    double_t m12;
    double_t m23;
};

