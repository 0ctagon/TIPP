
#include "TLorentzVector.h"

void test() {
// example of use of TGenPhaseSpace
//Author: Valerio Filippini

   if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");

   TLorentzVector b0(0.0, 0.0, 0.65, 1.538);

   RooRealVar phi2("phi2","angle between 2nu in the COM of B",50,0,5884);
   std::cout<<phi2.getMax()<<" "<<phi2.getMin()<<std::endl;
   //(Momentum, Energy units are Gev/C, GeV)
   Double_t masses[3] = { 0.938, 0., 0.} ;

   TGenPhaseSpace event;
   event.SetDecay(b0, 3, masses);

   TH2F *h2 = new TH2F("h2","h2", 50,1.1,1.8, 50,1.1,1.8);

   for (Int_t n=0;n<100000;n++) {
      Double_t weight = event.Generate();

      TLorentzVector *pProton = event.GetDecay(0);

      TLorentzVector *pPip    = event.GetDecay(1);
      TLorentzVector *pPim    = event.GetDecay(2);

      TLorentzVector pPPip = *pProton + *pPip;
      TLorentzVector pPPim = *pProton + *pPim;

      h2->Fill(pPPip.M2() ,pPPim.M2() ,weight);
   }
   h2->Draw();
}