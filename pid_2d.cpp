#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"

using namespace std;

int main(){
  int mc_steps=1e2;

  double amu=938;     //MeV/c^2
  double length = 50; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = .7; // truncation factor

  Track pion{140,500,length,segment,factor};
  Track p{amu,100,length,segment,factor};
  Track d{amu*2,100,length,segment,factor};
  Track t{amu*3,100,length,segment,factor};

  pion.SetInvXSec("P10M0invw_31623.inv");
  p.SetInvXSec("P10M0invw_31623.inv");
  d.SetInvXSec("P10M0invw_31623.inv");
  t.SetInvXSec("P10M0invw_31623.inv");

     cout << "Bgamma is " << pion.Getbg() << endl;
    cout << "Mean free path " << pion.GetMpath() <<endl;

    TH2D *pi_hist;
    TH2D *p_hist;
    TH2D *d_hist;
    TH2D *t_hist;

    pi_hist = pion.GraphMomRange(mc_steps,100,50,1500);
    p_hist =     p.GraphMomRange(mc_steps,150,100,3000);
    d_hist =     d.GraphMomRange(mc_steps,150,200,3000);
    t_hist =     t.GraphMomRange(mc_steps,150,300,3000);

    TFile *out = new TFile("pid.root","RECREATE");
    TH1D *f;
    TCanvas *c1 = new TCanvas(1);

    t_hist->Add(pi_hist);
    t_hist->Add(p_hist);
    t_hist->Add(d_hist);

  t_hist->Draw("colz");
  t_hist->RebinX(4);
  t_hist->GetYaxis()->SetRangeUser(0,16);
    t_hist->GetXaxis()->SetRangeUser(0,3000);

    t_hist->Write();
    pi_hist->Write();
    p_hist->Write();
    d_hist->Write();
    
   c1->SaveAs("pid.png");

  return 0;
}
