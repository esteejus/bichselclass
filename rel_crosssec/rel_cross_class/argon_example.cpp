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

  argon.GetRelCrossSection(3.16);
  TGraph *cross = argon.DrawCrossSection(true);
  cout<<"Moment 0 "<<argon.GetMoment(0,3.16)<<endl;;
  cout<<"Moment 1    "<<argon.GetMoment(1,3.16)/1e3<<endl;;
  cout<<"BETHE BLOCH "<<  argon.GetBetheBloch(3.16)<<endl;
  
  TGraph *bichsel = argon.DrawBichselSeg(938,2964.08,.5,1e3,0,5e3);
  //  TGraph *ruth = argon.DrawRutherford(1e2,10,500,true);

  TCanvas *c7 = new TCanvas("c7","c7",1);
  c7 -> SetLogx();
  bichsel -> Draw();
  c7 -> SaveAs("bichsel_seg.png");
  
  TCanvas *c6 = new TCanvas("c6","c6",1);
  //  c6 -> SetLogy();
  c6 -> SetLogx();
  // cross -> GetYaxis() -> SetRangeUser(1e-6,.08);
  //  cross -> GetXaxis() -> SetRangeUser(8,100);
  cross -> GetXaxis() -> SetLimits(8,100);
  //  ruth -> SetLineColor(2);
  //  ruth ->Draw("ALO");
  cross -> Draw("ALO");
  c6 -> SetLogx();
  c6 -> SaveAs("cross_section.png");

  return 0;
}
