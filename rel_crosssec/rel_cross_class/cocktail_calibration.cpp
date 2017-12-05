#include "dielectric.h"
#include "TFile.h"
#include <gsl/gsl_spline.h>

using namespace std;

int main()
{
  Dielectric argon{1.66e-3,38.,18,"argon"};
  argon.SetPhotoCross("../argon_10_500_Sakamoto_ext.dat");

  //Set real and imaginary dielectric tables for faster computaiton
  if(argon.SetRealTable("real_argon.dat") == false)
    {
      argon.GetRealDielectric();
      argon.WriteToFile("real");
    }
  if(argon.SetImgTable("img_argon.dat") == false)
    {
      argon.GetImgDielectric();
      argon.WriteToFile("img");
    }

   double amu     = 931.5;     //MeV/c^2 from Nuclear Wallet cards
   int    ns      = 108; // number of segments
   double segment = 1.2; // [cm] segment analyzed
   double factor  = .7; // truncation factor

  //masses from nuclear wallet cards
    double p_mass = amu+7.2889;
    double d_mass = amu*2+13.136;
    double t_mass = amu*3+14.95;
    double he3_mass = amu*3+14.931;
    double he4_mass = amu*4+2.425;
  
    double p_mom = 1700;
    double d_mom = 1700;
    double t_mom = 1700;
    //    double he3_mom = 3236.4;
    //    double he4_mom = 3226.35;

    TGraph *p_dist = argon.DrawBichselSeg(p_mass,p_mom,segment,1e4,0,1e6,0);
    TGraph *d_dist = argon.DrawBichselSeg(d_mass,d_mom,segment,1e4,0,1e6,0);
    TGraph *t_dist = argon.DrawBichselSeg(t_mass,t_mom,segment,1e4,0,1e6,0);
    
    TH1D *p_c_dist = argon.GetMCdist(p_dist,"p_c_dist",ns,factor,1e6);
    TH1D *d_c_dist = argon.GetMCdist(d_dist,"d_c_dist",ns,factor,1e6);
    TH1D *t_c_dist = argon.GetMCdist(t_dist,"t_c_dist",ns,factor,1e6);

    TCanvas *c1 = new TCanvas("c1","c1",1);
    d_c_dist->SetLineColor(2);
    t_c_dist->SetLineColor(4);

    p_c_dist ->Draw();
    d_c_dist->Draw("");
    t_c_dist->Draw("same");
    
    c1 -> SaveAs("c_dist_cocktail.png");

    TFile *out = new TFile("cocktail300.root","RECREATE");
       p_c_dist -> Write();
       d_c_dist -> Write();
       t_c_dist -> Write();
       out->Close();
    
  return 0;
}
