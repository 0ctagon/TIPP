
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

void GraphEffEmiss(int N)
{
    double_t Veff1[N],Vpur1[N],Vrej1[N];
    double_t Veff2[N],Vpur2[N],Vrej2[N];
    double_t Veff3[N],Vpur3[N],Vrej3[N];

    ifstream Vect("sigEmiss.txt");
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
    GraphPurEff1->GetHistogram()->SetMaximum(0.06);
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

    TCanvas* c = new TCanvas("comparisonsigEmiss","comparison sigEmiss",900,450) ;
	c->Divide(2,1);
	c->cd(1);
	GraphPurEff1->Draw("APL") ;
    GraphPurEff2->Draw("same,PL") ;
    GraphPurEff3->Draw("same,PL") ;
    c->cd(2);
	GraphRejEff1->Draw("APL") ;
    GraphRejEff2->Draw("same,PL") ;
    GraphRejEff3->Draw("same,PL") ;

    TLegend *leg = new TLegend( .58, .65, .9, .9, "Emiss Resolution");
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->AddEntry( GraphPurEff1, "5\%", "lp");
    leg->AddEntry( GraphPurEff2, "10\%", "lp");
    leg->AddEntry( GraphPurEff3, "30\%", "lp");
    c->cd(1);leg->Draw();
    Vect.close();
}