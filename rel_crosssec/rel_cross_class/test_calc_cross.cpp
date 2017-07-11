#include "dielectric.h"

using namespace std;

int main()
{
  Dielectric d1{1.66e-3,38.,18};
  d1.SetPhotoCross("../argon_10_500_Sakamoto.dat");
  //  d1.GetImgDielectric(1e4,0,1.e3);
  //  d1.GetRealDielectric(1e4,0,1.e3);

  //  d1.GetRelCrossSection(3.6,1e4,0,1e3);

  Dielectric d2{.6669e-3,16.043,10};
  d2.SetPhotoCross("../ch4.dat");
  //  d2.GetImgDielectric(1e4,0,1.e3);
  //  d2.GetRealDielectric(1e4,0,1.e3);

  Dielectric d3 = d1.MixGas(d1,d2,.9,.1,1e5,0,1e4);
  //  d3.GetImgDielectric(1e4,0,1.e3);
  //  d3.GetRealDielectric(1e4,0,1.e3);
  //  cout<<"z "<<d3.GetZ()<<" "<<d3.GetMolMass()<<" "<<d3.GetAtomcm3()<<endl;

    vector<double> energy = d3.GetEnergyVec();
    vector<double> value = d3.GetPhotoCrossVec();

    cout<<"size is "<<  energy.size()<<endl;
    //    for(int i =0;i<energy.size();i++)
      //    cout<<energy.at(i)<<" "<<value.at(i)<<endl;

  
    TGraph *photo = d3.DrawPhotoCross(1e4,0,1.e3);
    //    TGraph *photo = d1.DrawPhotoCross(1e5,0,1.e4);
  /*
  TGraph *ar_img = d1.DrawImaginary();
  TGraph *ar_real = d1.DrawReal();
  TGraph *cross = d1.DrawCrossSection(true);
  TGraph *ruth = d1.DrawRutherford(1e2,10,500,true);

  TGraph *mix_img = d3.DrawImaginary();
  */
  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1 -> SetLogy();
  c1 -> SetLogx();
  photo -> Draw();
  c1 -> SaveAs("testplotcross.png");
  /*
  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->SetLogy();
  c3->SetLogx();
  ar_img -> Draw();
  c3 -> SaveAs("ar_img.png");

  TCanvas *c5 = new TCanvas("c5","c5",1);
  c5 -> SetLogx();
  ar_real -> Draw();
  c5 -> SaveAs("ar_real.png");

  TCanvas *c6 = new TCanvas("c6","c6",1);
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
