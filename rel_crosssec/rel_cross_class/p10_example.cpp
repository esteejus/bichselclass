#include "dielectric.h"

using namespace std;

int main()
{
  Dielectric argon{1.66e-3,38.,18,"argon"};//Argon @ 1 atm
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

  Dielectric ch4{.667e-3,16.,10,"ch4"};//CH4 @ 1 atm
  ch4.SetPhotoCross("../ch4.dat");

  //Set real and imaginary dielectric tables for faster computaiton
  if(ch4.SetRealTable("real_ch4.dat") == false)
    {
      ch4.GetRealDielectric();
      ch4.WriteToFile("real");
    }
  if(ch4.SetImgTable("img_ch4.dat") == false)
    {
      ch4.GetImgDielectric();
      ch4.WriteToFile("img");
    }  

Dielectric p10 = p10.MixGas(argon,ch4,.9,.1);
 p10.SetName("p10");
 cout<<"Density of "<<p10.GetName()<<" "<<p10.GetDensity()<<" Mol mass "<<p10.GetMolMass()<<" Z "<<p10.GetZ()<<endl;
  
  //Set real and imaginary dielectric tables for faster computaiton
  if(p10.SetRealTable("real_p10.dat") == false)
    {
      p10.GetRealDielectric();
      p10.WriteToFile("real");
    }
  if(p10.SetImgTable("img_p10.dat") == false)
    {
      p10.GetImgDielectric();
      p10.WriteToFile("img");
    }  

  TGraph *ch4_real = ch4.DrawReal();
  TGraph *argon_real = argon.DrawReal();
  TGraph *ch4_img = ch4.DrawImaginary();
  TGraph *argon_img = argon.DrawImaginary();
  TGraph *p10_real = p10.DrawReal();
  TGraph *p10_img = p10.DrawImaginary();
  p10_real -> SetLineColor(1);
  p10_img -> SetLineColor(1);
  argon_real -> SetLineColor(4);
  argon_img -> SetLineColor(4);
  ch4_real -> SetLineColor(2);
  ch4_img -> SetLineColor(2);

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
  
  //  double bgamma = .4;
  double p_mom = 935;
  double d_mom = 1700;
  double t_mom = 1700;
  double he3_mom = 1700*2;
  double he4_mom = 1700*2;
  
  //  TGraph *  bichsel_argon = argon.DrawBichselSeg(p_mass,p_mom,segment,npoints,0,14000e3);
  //TH1D *  c_dist_argon = argon.GetMCdist(bichsel_argon,"c_dist_argon",segment,ns,factor,1e6);
 
TGraph *  bichsel_p10 = p10.DrawBichselSeg(t_mass,t_mom,segment,npoints,0,14000e3);
 TH1D *  c_dist_p10 = p10.GetMCdist(bichsel_p10,"c_dist_p10",segment,ns,factor,1e6); 
 
 /*
 TCanvas *c2 = new TCanvas("c2","c2",1);
 c2->SetLogx();
 // bichsel_argon->GetXaxis()->SetLimits(1e3,1e4);
 // bichsel_argon->Draw();
 // bichsel_p10->GetXaxis()->SetLimits(1e3,6e3);
 // bichsel_p10->SetLineColor(2);
 // bichsel_p10->Draw("ALO");
 c_dist_p10->Draw("ALO");
 c2->SaveAs("p10_argon_bichsel.png");

 cout<<"mean value "<<c_dist_p10->GetMean()<<endl;
 /*TCanvas *c7 = new TCanvas("c7","c7",1);
  c7 -> SetLogx();
  //  argon_real->Draw("ALO");
  //  p10_real->Draw("same LO");
  //  ch4_real->Draw("same LO");
  //  c7 -> SaveAs("real_gas.png");

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1 -> SetLogx();
  c1 -> SetLogy();
  //  ch4_img -> GetXaxis()-> SetRangeUser(.1,1e4);
  //  ch4_img->Draw("ALO");
//  p10_img ->Draw("same LO");
  // argon_img->Draw("same LO");

  c1 -> SaveAs("img_gas.png");
 */
 ofstream outfile;
  outfile.open("bichsel_theory.dat");
  
  for(int i = 1; i < bichsel_p10->GetN();i++)
    {
      double mean = 0.;
      double value =  0.;
      bichsel_p10->GetPoint(i,mean,value);
      if(value > 1e-10 && value > 0)
	outfile<<mean/1000<<"\t"<<value<<endl;
    }

  ofstream cfile;
  cfile.open("cstrag_theory.dat");
  
  for(int i = 1; i <= c_dist_p10->GetXaxis()->GetNbins();i++)
    {
      double mean = c_dist_p10->GetBinCenter(i);
      double value = c_dist_p10->GetBinContent(i);

      if(value > 1e-10 && value > 0)
	cfile<<mean/1000<<"\t"<<value<<endl;
    }
  



  return 0;
}
