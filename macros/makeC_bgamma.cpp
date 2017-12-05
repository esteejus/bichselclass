#include "track.h"
#include "bichsel.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include <fstream>

using namespace std;

int main(){
  int mc_steps=1e5;

  double amu=931.5;     //MeV/c^2 from Nuclear Wallet cards
  double length = 100; // [cm] length of track
  double segment = 1.2; // [cm] segment analyzed
  double factor = .7; // truncation facto

  int npoints = 1;
  //masses from nuclear wallet cards
    double mass = amu+7.289;//mass of proton
    ofstream out;
    out.open("C_avg_6.dat");

     double bgamma = .017;
     double mom = bgamma*mass;
     Track part{mass,mom,length,segment,factor};
     part.SetInvXSec("P10M0invw_31623.inv");
     part.SetM0Table("m0.dat");
     part.GetCArray(1e4);
     out<<bgamma<<"\t"<<part.GetCavg()<<"\t"<< part.GetCsigma()<<endl;

  return 0;
}



