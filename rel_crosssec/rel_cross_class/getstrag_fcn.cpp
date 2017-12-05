#include "dielectric.h"

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
  
  int npoints = 1e3;
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
  
  double p_mom = 925;
  double d_mom = 1700;
  double t_mom = 1700;
  double he3_mom = 1700*2;
  double he4_mom = 1700*2;
  
 TGraph *  bichsel = argon.DrawBichselSeg(d_mass,d_mom,segment,npoints,0,14000e3);
 TH1D *  c_dist = argon.GetMCdist(bichsel,"c_dist",segment,ns,factor,1e6);  //  TGraph 

  TCanvas *c7 = new TCanvas("c7","c7",1);
  c7 -> SetLogx();
  bichsel -> Draw();
  c7 -> SaveAs("bichsel_seg.png");
  

  ofstream outfile;
  outfile.open("bichsel_theory.dat");
  
  /* for(int i = 1; i < bichsel->GetN();i++)
    {
      double mean = 0.;
      double value =  0.;
      bichsel->GetPoint(i,mean,value);
      if(value > 1e-10 && value > 0)
	outfile<<mean/1000<<"\t"<<value<<endl;
    }
  */
  ofstream cfile;
  cfile.open("cstrag_theory.dat");
  
  for(int i = 1; i <= c_dist->GetXaxis()->GetNbins();i++)
    {
      double mean = c_dist->GetBinCenter(i);
      double value = c_dist->GetBinContent(i);

      if(value > 1e-10 && value > 0)
	cfile<<mean/1000<<"\t"<<value<<endl;
    }
  
  
 return 0;
}
