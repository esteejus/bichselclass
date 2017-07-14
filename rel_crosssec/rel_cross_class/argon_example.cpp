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

  argon.GetRelCrossSection(.316);
  TGraph *cross = argon.DrawCrossSection(true);
  cout<<"Moment 0 "<<argon.GetMoment(0)<<endl;;
  cout<<"Moment 1    "<<argon.GetMoment(1)/1e3<<endl;;
  cout<<"BETHE BLOCH "<<  argon.GetBetheBloch(.316)<<endl;
  
    //  TGraph *ruth = argon.DrawRutherford(1e2,10,500,true);

  TCanvas *c6 = new TCanvas("c6","c6",1);
  c6 -> SetLogx();
  c6 -> SetLogy();
  //  ruth -> GetYaxis() -> SetRangeUser(1e-6,2);
  //  ruth -> SetLineColor(2);
  //  ruth ->Draw("ALO");
  cross -> Draw("ALO");
  c6 -> SaveAs("cross_section.png");

  return 0;
}
