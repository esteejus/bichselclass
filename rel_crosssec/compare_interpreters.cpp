#include <gsl/gsl_integration.h>
//#include <gsl/gsl_erron.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "spline.h"
#include <cmath>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
using namespace std;

double photo_cross_interp(double, void *);

std::vector<double> photoenergy; //array to store values of cumulative dist
std::vector<double> photovalue; //array to store values of cumulative dist

gsl_interp_accel *acc = gsl_interp_accel_alloc();
gsl_spline *photo_cross_table_linear;
gsl_spline *photo_cross_table_cspline;
gsl_spline *photo_cross_table_akima;

void SetPhotoCross (const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      photoenergy.push_back(d_energy);
      photovalue.push_back(d_value);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  double *x_a = photoenergy.data();
  double *y_a = photovalue.data();
  int size_a = photoenergy.size();
  photo_cross_table_linear = gsl_spline_alloc(gsl_interp_linear,size_a);
  photo_cross_table_cspline = gsl_spline_alloc(gsl_interp_cspline,size_a);
  photo_cross_table_akima = gsl_spline_alloc(gsl_interp_akima,size_a);

  gsl_spline_init(photo_cross_table_linear,x_a,y_a,size_a);
  gsl_spline_init(photo_cross_table_cspline,x_a,y_a,size_a);
  gsl_spline_init(photo_cross_table_akima,x_a,y_a,size_a);
  
  return;
}


int main(){

  ofstream output;
  SetPhotoCross("argon_west.dat");

  int npoints = 1e6;

  TGraph interp_linear = TGraph(npoints);
  TGraph interp_cspline = TGraph(npoints);
  TGraph interp_akima = TGraph(npoints);
  
  double emin = 16.;
  double emax = 4e3;
  double energystep = (emax-emin)/npoints;

  int num = photoenergy.size();
  TGraph cross = TGraph(num);
  
  for(int i=0;i<photoenergy.size();++i){
    double x = photoenergy.at(i);
    double y = photovalue.at(i);
    cross.SetPoint(i,x,y);
  }

  for(int i=0;i<npoints;++i){
    double x = emin + energystep*i;

    double lin_v    = gsl_spline_eval(photo_cross_table_linear,x,acc);
    double spline_v = gsl_spline_eval(photo_cross_table_cspline,x,acc);
    double akima_v  = gsl_spline_eval(photo_cross_table_akima,x,acc);
 
    interp_linear.SetPoint(i,x,lin_v);
    interp_cspline.SetPoint(i,x,spline_v);
    interp_akima.SetPoint(i,x,akima_v);
  }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->SetLogx();
  c1->SetLogy();
  cross.GetXaxis()->SetRangeUser(3.35e3,3.45e3);
  cross.SetMarkerSize(4);
  cross.Draw();
  
  interp_linear.SetMarkerColor(1);
    interp_linear.Draw("same PO");

  interp_cspline.SetMarkerColor(2);
    interp_cspline.Draw("same PO");

  interp_akima.SetMarkerColor(4);
    interp_akima.Draw("same PO");

  c1->SaveAs("interpcompare.png");

  
  return 0;
  
}



