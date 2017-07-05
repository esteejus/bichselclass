#include "dielectric.h"

using namespace std;

int main()
{
  Dielectric d1{1.66e-3,38.};
  d1.SetPhotoCross("../argon_10_500_Sakamoto.dat");
  d1.GetImgDielectric(1e4,0,1.e3);
  d1.GetRealDielectric(2e3,0,1.e3);
// d1.SetPhotoCross("../ch4.dat");
  TGraph *photo = d1.DrawPhotoCross(1e5,0,1.e4);
  TGraph *ar_img = d1.DrawImaginary();
  TGraph *ar_real = d1.DrawReal();

  Dielectric d2{1.66e-3,38.};
  d2.SetPhotoCross("../ch4.dat");
  d2.GetImgDielectric(1e4,0,1.e3);
  TGraph *photo_ch4 = d2.DrawPhotoCross(1e5,0,1.e4);
  TGraph *ch4_img = d2.DrawImaginary();

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1 -> SetLogy();
  c1 -> SetLogx();
  photo -> Draw();
  c1 -> SaveAs("testplotcross.png");

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2 -> SetLogy();
  c2 -> SetLogx();
  photo_ch4 -> GetYaxis()-> SetRangeUser(1e-6,50);
  photo_ch4 -> Draw();
  c2 -> SaveAs("testplotcross_ch4.png");

  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->SetLogy();
  c3->SetLogx();
  ar_img -> Draw();
  c3 -> SaveAs("ar_img.png");

  TCanvas *c4 = new TCanvas("c4","c4",1);
  c4->SetLogy();
  c4->SetLogx();
  ch4_img -> Draw();
  c4 -> SaveAs("ch4_img.png");

  TCanvas *c5 = new TCanvas("c5","c5",1);
  c5 -> SetLogx();
  ar_real -> Draw();
  c5 -> SaveAs("real.png");

  return 0;
}
