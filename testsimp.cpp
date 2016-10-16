#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"
#include "TGraph.h"

using namespace std;

int main(){
  int mc_steps=1e3;

  double amu=938;     //MeV/c^2
  double length = 50; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = .7; // truncation factor

  Track pion{140,50,length,segment,factor};
  pion.SetInvXSec("P10M0invw_31623.inv");

  TH1 *f;
  f = pion.Drawfdist(1e4,20);
  //  pion.HistArray(.1,0,20);
  
    std::vector<std::vector<double> > hist = pion.SimpsonNInt(.02,0,20);
    int num = hist[0].size();
    TGraph *g = new TGraph(num);
  
  for(int i=0;i<hist[0].size();++i){
           g->SetPoint(i,hist[0][i],hist[1][i]);
  }
  
  TCanvas *c1 = new TCanvas(1);
  g->Draw();

  c1->SaveAs("cumulative.jpg");

  return 0;
}
