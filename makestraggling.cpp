#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"
#include "TFile.h"
#include <fstream>

using namespace std;

int main(){
  int mc_steps=1e5;

  double amu=931.5;     //MeV/c^2 from Nuclear Wallet cards
  double length = 100; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = 1.; // truncation facto

    string name = "t";

  //masses from nuclear wallet cards
  double p_mass = amu+7.289;
  double d_mass = amu*2+13.136;
  double t_mass = amu*3+14.95;
  double he3_mass = amu*3+14.931;
  double he4_mass = amu*4+2.425;

  //Uncommentfor 300 MeV/A beam  
  double d_mom = 1621.05;
  double t_mom = 1612.35;
  double he3_mom = 3236.4;
  double he4_mom = 3226.35;

  double mass = -99.;
  double mom  = -999.;
  
  cout<<"Enter particle type (p,d,t,he3,he4):"<<endl;
  cin>>name;
  cout<<"Enter momentum [MeV/c]:"<<endl;
  cin>>mom;
 
    if(name=="d")        mass = d_mass;
    else if(name=="p")   mass = p_mass;
    else if(name=="t")   mass = t_mass;
    else if(name=="he3") mass = he3_mass;
    else if(name=="he4") mass = he4_mass;
    else cout<<"does not match known particles"<<endl;


    cout<<"particle mass is "<<mass<<endl;
  
  Track part{mass,mom,length,segment,factor};
  part.SetInvXSec("P10M0invw_31623.inv");
  TH1D *dedx_hist;

 dedx_hist = part.Drawfdist(mc_steps,20);
 cout<<dedx_hist->Integral();
 // dedx_hist->Scale(1./dedx_hist->Integral());

  TCanvas *c1 = new TCanvas(1);
  dedx_hist->Draw();
  c1->SetLogy();
  TString c_name = name + "_dedx_fdist.png";
  c1->SaveAs(c_name);

  ofstream out;
  out.open(name+"_dedx.dat");
  out<<"Particle "<<name<<" momenum [MeV/c] "<< mom << " mass [MeV/c^2] "<<mass;
  int num_binsx = dedx_hist->GetNbinsX();
  for(int i = 1; i<= num_binsx; i++){
    out<<dedx_hist->GetBinCenter(i)<<"\t"<<dedx_hist->GetBinContent(i)<<endl;  
  }
  out.close();  
 
  TString r_name = name + "dedx.root";
  TFile *f = new TFile(r_name,"RECREATE");
  //  dedx_hist->Scale(1./dedx_hist->Integral("width"));
  dedx_hist->Write();
  
  
  return 0;
}
