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

  int nsteps = 25;
  double bgamma_l = 1.1;
  double bgamma_h = 4.;
  double bgamma_log_step = (log(bgamma_h)-log(bgamma_l))/nsteps;

  ofstream outfile;
  outfile.open("bgammadEdx.dat");
  
  for( int i = 0 ;i < nsteps; i++)
    {
      double bgamma = bgamma_log_step*i + log(bgamma_l);
      bgamma = exp(bgamma);
      cout<<"Beta gamma is "<<bgamma<<endl;      
      argon.GetRelCrossSection(bgamma);
  double dummy_mass = 938;
  double mom = bgamma*dummy_mass;
  TGraph *  bichsel = argon.DrawBichselSeg(dummy_mass,mom,1.2,1e3,0,14000e3);
  TH1D *  c_dist = argon.GetMCdist(bichsel,"c_dist",1.2,90,.7,1e6);
  outfile<<bgamma<<"\t"<<  c_dist -> GetMean()<<endl;

    }



  return 0;
}
