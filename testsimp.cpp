#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"

using namespace std;

int main(){
  int mc_steps=1e3;

  double amu=938;     //MeV/c^2
  double length = 50; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = .7; // truncation factor

  Track pion{140,100,length,segment,factor};
  pion.SetInvXSec("P10M0invw_31623.inv");

  TH1 *f;
  f = pion.Drawfdist(1,20);
  pion.HistArray(.1,0,20);
  
  return 0;
}
