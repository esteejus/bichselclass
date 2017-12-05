#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace std;

int main(){
  int mc_steps=5e3;

  double amu=938;     //MeV/c^2
  double length = 50; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = .7; // truncation factor

  Track pion{140,500,length,segment,factor};
  Track p{amu,100,length,segment,factor};
  Track d{amu*2,100,length,segment,factor};
  Track t{amu*3,100,length,segment,factor};

  TFile *out = new TFile("pid_high.root","RECREATE");
  
  pion.SetInvXSec("P10M0invw_31623.inv");

     cout << "Bgamma is " << pion.Getbg() << endl;
    cout << "Mean free path " << pion.GetMpath() <<endl;

    TH2D *pi_hist_high;
    
    pi_hist_high = pion.GraphMomRange(mc_steps,400,200,2000);
    pi_hist_high->Write();

  TCanvas *c1 = new TCanvas(1);
    pi_hist_high->Draw("colz");
    //    pi_hist_high->RebinX(4);
    pi_hist_high->GetYaxis()->SetRangeUser(0,16);
    pi_hist_high->GetXaxis()->SetRangeUser(0,2500);

   c1->Draw();
   c1->SaveAs("pid_hi.png");
   c1->Write();

   out->Close();
   
  return 0;
}
