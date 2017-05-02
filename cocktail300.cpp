#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace std;

int main(){
  int mc_steps=1e3;

  double amu=931.5;     //MeV/c^2 from Nuclear Wallet cards
  double length = 1188; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = 1.; // truncation factor

  //masses from nuclear wallet cards
    double d_mass = amu*2+13.136;
  double t_mass = amu*3+14.95;
  double he3_mass = amu*3+14.931;
  double he4_mass = amu*4+2.425;
  
  double d_mom = 1621.05;
  double t_mom = 1612.35;
  double he3_mom = 3236.4;
  double he4_mom = 3226.35;
    
  Track d{d_mass,d_mom,length,segment,factor};
  Track t{t_mass,t_mom,length,segment,factor};
  Track he3{he3_mass,he3_mom,length,segment,factor};
  Track he4{he4_mass,he4_mom,length,segment,factor};

  d.SetInvXSec("P10M0invw_31623.inv");
  t.SetInvXSec("P10M0invw_31623.inv");
  he3.SetInvXSec("P10M0invw_31623.inv");
  he4.SetInvXSec("P10M0invw_31623.inv");

  TFile *out = new TFile("cocktail300.root","RECREATE");
  TH1D *f[4];

  f[0] = d.Drawfdist(mc_steps,20);
  f[1] = t.Drawfdist(mc_steps,20);
  f[2] = he3.Drawfdist(mc_steps,20);
  f[3] = he4.Drawfdist(mc_steps,20);

  f[0]->SetName("d");
  f[1]->SetName("t");
  f[2]->SetName("he3");
  f[3]->SetName("he4");
  
  TCanvas *c1 = new TCanvas(1);
  f[0]->SetLineColor(1);
  f[1]->SetLineColor(2);
  f[2]->SetLineColor(3);
  f[3]->SetLineColor(4);

  f[0]->Draw();
  f[1]->Draw("same");
  f[2]->Draw("same");
  f[3]->Draw("same");

  f[0]->Write();
  f[1]->Write();
  f[2]->Write();
  f[3]->Write();
  
   c1->SaveAs("cocktail300.png");

  return 0;
}
