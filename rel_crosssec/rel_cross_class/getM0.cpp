#include "dielectric.h"
#include <fstream>

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
  int npoints = 1e2;
  double bgamma_l = .01;
  double bgamma_h = 3.5;
  double beta_step = (bgamma_h-bgamma_l)/npoints;

  double log_step = (log(bgamma_h)-log(bgamma_l))/npoints;

  ofstream output;
  output.open("m0.dat");
  
  for(int i = 0; i <= npoints; i++){
    double ep = log_step*i+log(bgamma_l);
    double bgamma = exp(ep);
    //    double bgamma = beta_step*i + bgamma_l;
    argon.GetRelCrossSection(bgamma);
    double value = argon.GetMoment(0,bgamma);
    output<<bgamma<<"\t"<<value<<endl;;
    cout<<"Beta gamma "<<bgamma<<" M0 "<<value<<endl;
  }

  return 0;
}
