#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"

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

