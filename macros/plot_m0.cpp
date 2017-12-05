#include "bichsel.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;


int main(){

  Bichsel b1(1);
  b1.SetM0Table("m0.dat");

  int npoints = 1e3;
  double bgamma_l = .3;
  double bgamma_h = 8;

  double bgamma_step = (bgamma_h-bgamma_l)/npoints;
  TGraph * bgamma_graph = new TGraph(npoints);
  for(int i = 0;i < npoints; i++)
    {
      double bgamma = bgamma_step*i + bgamma_l;
      double value = b1.InterpolateM0(bgamma);
      bgamma_graph -> SetPoint(i,bgamma,value);
    }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  bgamma_graph -> Draw();
  c1 -> SaveAs("m0_graph.png");
    
  double mass,mom,bgamma;
  cout<<"Enter mass "<<endl;
  cin>>mass;
  cout<<"Enter mom"<<endl;
  cin>>mom;
  bgamma = mom/mass;
  double value = b1.InterpolateM0(bgamma);
  cout<<"dE/dx is "<<endl;
  cout<<value<<endl;
  
  return 0;
}
