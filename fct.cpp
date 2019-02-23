#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include <iostream>
using namespace std;

//Define the different variables inside a Particle object
struct Particle{
   	double_t E;
    double_t P;
    double_t Px;
    double_t Py;
    double_t Pz;
    double_t Costheta;
    double_t Phi;
};

//Difine the different variables used to get the Dalitz plot
struct Dalitz{
    double_t m12;
    double_t m23;
};

//Return the boost vector given a 4-vector
TVector3 Boost(TLorentzVector A)
{
    TVector3 beta(A.Px()/A.E(),A.Py()/A.E(),A.Pz()/A.E());
    return beta;
}

//Return a 4-vector in the lab frame
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

//Fill list and histogram with efficiency values
void Efficiency(double_t* Veff, TH1* THv, TH1* hist, int Nbin, double_t xmin, double_t xmax)
{
    double_t full_integral = hist->Integral(), cut_integral;

    for(int i=0;i<Nbin+1;i+=1)
    {
        cut_integral = hist->Integral(i,Nbin)/full_integral;
        Veff[i]=cut_integral;
        THv->SetBinContent(i,Veff[i]);
    }
}

//Fill list and histogram with rejection values
void Rejection(double_t* Vrej, TH1* THv, TH1* hist, int Nbin, double_t xmin, double_t xmax)
{
    double_t cut_integral;
    double_t full_integral = hist->Integral();

    for(int i=0;i<Nbin+1;i+=1)
    {
        cut_integral = hist->Integral(i,Nbin)/full_integral;
        Vrej[i]=1-cut_integral;
        THv->SetBinContent(i,Vrej[i]);
    }
}

//Fill list and histogram with purity values
void Purity(double_t* Vpur, TH1* THv, TH1* histbkg, TH1* histE, int Nbin, double_t xmin, double_t xmax)
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
        THv->SetBinContent(i,Vpur[i]);
    }
}

//Return the maximum significance value or the cut value corresponding to this maximum and fill a histogram (TH1sign) with all the significance values
double_t SignificanceTH(TH1* TH1sign, TH1* histEff, TH1* histRej, int Nbin, int Nbkg, int f, double_t xmin, double_t xmax, bool cut)
{
    double_t Vsignificance[Nbin];

    double_t eff,rej;

    for(int i=0;i<Nbin+1;i+=1)
    {
        eff = histEff->GetBinContent(i-1);
        rej = 1-histRej->GetBinContent(i-1);

        if((sqrt(eff+f*rej))<0.001)
        {
            Vsignificance[i]=0.0;
        }
        else
        {
            Vsignificance[i] = (sqrt(Nbkg/f)*eff)/(sqrt(eff+f*rej));
        }

        TH1sign->SetBinContent(i,Vsignificance[i]);
    }

    double_t value;
    if(cut)
    {
        value = ((xmax-xmin)/Nbin)*(TH1sign->GetMaximumBin()-0.5);
    }
    else
    {
        value = TH1sign->GetMaximum();
    }
    
    return value;
}


TH1D *TH1sign = new TH1D("TH1sign","purity",100,0.0,7.0);

//Return the maximum significance value or the cut value corresponding to this maximum
double_t Significance(TH1* histEff, TH1* histRej, int Nbin, int Nbkg, int f, double_t xmin, double_t xmax, bool cut)
{
    TH1sign->BufferEmpty();
    double_t Vsignificance[Nbin];

    double_t eff,rej;

    for(int i=0;i<Nbin+1;i+=1)
    {
        eff = histEff->GetBinContent(i-1);
        rej = 1-histRej->GetBinContent(i-1);

        if((sqrt(eff+f*rej))<0.001)
        {
            Vsignificance[i]=0.0;
        }
        else
        {
            Vsignificance[i] = (sqrt(Nbkg/f)*eff)/(sqrt(eff+f*rej));
        }

        TH1sign->SetBinContent(i,Vsignificance[i]);
    }

    double_t value;
    if(cut)
    {
        value = ((xmax-xmin)/Nbin)*(TH1sign->GetMaximumBin()-0.5);
    }
    else
    {
        value = TH1sign->GetMaximum();
    }

    return value;
}

