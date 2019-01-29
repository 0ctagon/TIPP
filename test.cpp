
#include "TLorentzVector.h"
RooAbsPdf * fct(RooRealVar& x)
{
	RooGenericPdf * testt = new RooGenericPdf("fvgu","x*sqrt((pow(5.27931,2)-pow(x+0.493677,2))*(pow(5.27931,2)-pow(x-0.493677,2)))",x);
   //RooRealVar * csbkgmean = new RooRealVar("csbkgmean","CS bkg mean",0) ;
	//RooRealVar * csbkgsigma = new RooRealVar("csbkgsigma","CS bkg sigma ",1) ;
   //RooGaussian * testt = new RooGaussian("testt","gaussian PDF",x,*csbkgmean,*csbkgsigma) ;
   return testt;
}



void test() {
   RooRealVar x("x","x",0,0,4.7856330);
   RooAbsPdf * f = fct(x);
   RooAbsReal * g = f->createCdf(x);

   //m12.setVal(m12distribution(m12)->createCdf(m12)->asTF(m12)->GetX(m12));

   //std::cout<<g->getVal(0.5)<<std::endl;

   RooPlot* frame = x.frame() ;
	//frame->SetNdivisions(505,"X");
	RooPlot* frame2 = x.frame() ;
   
	f->plotOn(frame) ;
   g->plotOn(frame2) ;

   TF1 * gcumuhisto = g->asTF(x);

   int N=10000;
   double_t X[N],Y[N];
   double h=1.0/N;
   int j=0;

   for(double i=0;i<1;i+=h)
   {
      X[j]=gcumuhisto->GetX(i);
      //std::cout<<X[j]<<" "<<Y[j]<<" "<<h<<endl;
      Y[j]=i;
      j+=1;
   }
   TGraph *gr1 = new TGraph (N,Y,X);

   TCanvas* c2 = new TCanvas("bphysicsbkg","bphysics bkg",900,900) ;
	c2->Divide(2,2);
	c2->cd(1);
	gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
	c2->cd(2);
	gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
	c2->cd(3);
	gPad->SetLeftMargin(0.15) ; gr1->Draw() ;
   c2->cd(4);
	gPad->SetLeftMargin(0.15) ; //gcumuhisto->Draw() ;

   TH1D *m23c = new TH1D("m23c","Missing energy",Nbin,0,25);
    tree->Draw("Dalitz.m23_2>>m23c","","goff");
    TH1D *m12c = new TH1D("m12c","Missing energy",Nbin,0,30);
    tree->Draw("Dalitz.m12_2>>m12c","","goff");
    TH2D *testf = new TH2D("testf","idfvf",Nbin,0,25,Nbin,0,30);

    double_t x1=0,x2=0,k1,k2;
    double_t h1=25.0/Nbin,h2=30.0/Nbin;

    for(int i=0;i<Nbin;i++)
    {
        k1 = m23c->FindFixBin(x1);
        x2=0;
        for(int j=0;j<Nbin;j++)
        {
            k2 = m12c->FindFixBin(x2);
            testf->SetBinContent(k1,k2,1.0);
            x2+=h2;
        }
        x1+=h1;
    }
}