#include "dielectric.h"

using namespace std;

int main()
{
  Dielectric argon{1.66e-3,38.,18,"argon"};
  cout<<"BETHE BLOCH"<<  argon.GetBetheBloch(3.6)<<endl;
  cout<<"BETHE BLOCH"<<  argon.GetBetheBloch(.36)<<endl;
  argon.SetPhotoCross("../argon_10_500_Sakamoto.dat");
  argon.GetImgDielectric(1e4,0,1.e3);
  argon.GetRealDielectric(1e3,0,1.e3);
  argon.WriteToFile("real");
  argon.WriteToFile("img");

    
  //    argon.GetRelCrossSection(3.6,1e4,0,1e3);

  Dielectric ch4{.6669e-3,16.043,10,"ch4"};
  ch4.SetPhotoCross("../ch4.dat");
  ch4.GetImgDielectric(1e4,0,1.e3);
  ch4.GetRealDielectric(1e3,0,1.e3);
  ch4.WriteToFile("real");
  ch4.WriteToFile("img");

  Dielectric p10 = argon.MixGas(argon,ch4,.9,.1);
  p10.GetImgDielectric(1e4,0,1.e3);
  p10.GetRealDielectric(1e3,0,1.e3);
  p10.WriteToFile("real");
  p10.WriteToFile("img");
  p10.WriteToFile("photocross");
  
  p10.GetRelCrossSection(3.6,1e4,0,1e3);
  cout<<"Moment 0"<<p10.GetMoment(0)<<endl;;
  cout<<"Moment 0"<<p10.GetMoment(1)<<endl;;
  
  TGraph *photo_p10   = p10.DrawPhotoCross(1e5,0,1.e4);
  TGraph *photo_argon = argon.DrawPhotoCross(1e5,0,1.e4);
  TGraph *photo_ch4   = ch4.DrawPhotoCross(1e5,0,1.e4);
  

  TGraph *img_p10   = p10.DrawImaginary();
  TGraph *img_argon = argon.DrawImaginary();
  TGraph *img_ch4   = ch4.DrawImaginary();

  TGraph *real_p10   = p10.DrawReal();
  TGraph *real_argon = argon.DrawReal();
  TGraph *real_ch4   = ch4.DrawReal();

  //  TGraph *cross = argon.DrawCrossSection(true);
  //  TGraph *ruth = argon.DrawRutherford(1e2,10,500,true);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1 -> SetLogy();
  c1 -> SetLogx();
  photo_p10   -> SetLineColor(1);
  photo_argon -> SetLineColor(2);
  photo_ch4   -> SetLineColor(4);

  photo_p10   -> Draw("ALO");
  photo_argon -> Draw("LO");
  photo_ch4   -> Draw("LO");
  c1 -> SaveAs("allcross_sections.png");
  

  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->SetLogy();
  c3->SetLogx();

  img_p10   -> SetLineColor(1);
  img_argon -> SetLineColor(2);
  img_ch4   -> SetLineColor(4);

  img_p10   -> Draw("ALO");
  img_argon -> Draw("LO");
  img_ch4   -> Draw("LO");

  c3 -> SaveAs("all_img.png");

  TCanvas *c5 = new TCanvas("c5","c5",1);
  c5 -> SetLogx();
    real_p10   -> SetLineColor(1);
  real_argon -> SetLineColor(2);
  real_ch4   -> SetLineColor(4);

  real_p10   -> Draw("ALO");
  //  real_argon -> GetYaxis() -> SetRangeUser(-.0005,.0016);
  real_argon -> Draw("LO");
  real_ch4   -> Draw("LO");
  c5 -> SaveAs("all_real.png");

  /*  TCanvas *c6 = new TCanvas("c6","c6",1);
  c6 -> SetLogx();
  c6 -> SetLogy();
  //  ruth -> GetYaxis() -> SetRangeUser(1e-6,2);
  ruth -> SetLineColor(2);
  ruth ->Draw("ALO");
  cross -> Draw("same LO");
  c6 -> SaveAs("cross_section.png");

  TCanvas *c7 = new TCanvas("c7","c7",1);
  c7 -> SetLogx();
  mix_img -> Draw();
  c7 -> SaveAs("mixgas_img.png");
  */
  return 0;
}
