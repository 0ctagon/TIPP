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

void GraphEff(int N)
{
    double_t Veff1[N],Vpur1[N],Vrej1[N];
    double_t Veff2[N],Vpur2[N],Vrej2[N];
    double_t Veff3[N],Vpur3[N],Vrej3[N];

    ifstream Vect("vectors.txt");
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

    TGraph *GraphPurEff1 = new TGraph (N,Veff1,Vpur1);
    GraphPurEff1->SetTitle("Purity / Efficiency");
    GraphPurEff1->GetHistogram()->GetXaxis()->SetTitle("Efficiency");
    GraphPurEff1->GetHistogram()->GetYaxis()->SetTitle("Purity");
    TGraph *GraphRejEff1 = new TGraph (N,Veff1,Vrej1);
    GraphRejEff1->SetTitle("Rejection / Efficiency");
    GraphRejEff1->GetHistogram()->GetXaxis()->SetTitle("Efficiency");
    GraphRejEff1->GetHistogram()->GetYaxis()->SetTitle("Rejection");

    TGraph *GraphPurEff2 = new TGraph (N,Veff2,Vpur2);
    TGraph *GraphRejEff2 = new TGraph (N,Veff2,Vrej2);
    TGraph *GraphPurEff3 = new TGraph (N,Veff3,Vpur3);
    TGraph *GraphRejEff3 = new TGraph (N,Veff3,Vrej3);

    TCanvas* c = new TCanvas("comparisonsigEmiss","comparison sigEmiss",900,450) ;
	c->Divide(2,1);
	c->cd(1);
	GraphPurEff1->Draw("AL") ;
    GraphPurEff2->Draw("same") ;
    GraphPurEff3->Draw("same") ;
    c->cd(2);
	GraphRejEff1->Draw("AL") ;
    GraphRejEff2->Draw("same") ;
    GraphRejEff3->Draw("same") ;
    Vect.close();
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

