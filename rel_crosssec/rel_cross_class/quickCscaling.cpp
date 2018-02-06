#include "dielectric.h"
#include "TFile.h"
#include "TF1.h"
using namespace std;


int main()
{
  Dielectric argon{1.66e-3,38.,18,"argon"};
  //   TFile *f = new TFile("cstrag_d_1700_out.root");
    TFile *f = new TFile("cstrag_p_935_out.root");
  TH1D *theory = (TH1D*)f->Get("theory");
  TH1D *detreponse = (TH1D*)f->Get("detreponse");

  TGraph *scaling = argon.GetCScaling(theory,detreponse);
  TF1 *line = new TF1("line","[0]*x + [1]",0,100000);
  line->SetParameters(.1,1.4);
  TCanvas *c1 = new TCanvas("c1","c1",1);
  scaling ->Draw("APO");
  scaling->Fit("line");
  c1 -> SaveAs("cscale_theory_detresponse.png");
  
  return 0;
}
