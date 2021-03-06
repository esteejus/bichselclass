#include <gsl/gsl_integration.h>
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

double im_epsilon(double, void *);
double photo_cross(double, void *);
double photo_cross_interp(double, void *);
double LAP_photo_cross(double , void *);

std::vector<double> photoenergy; //array to store values of cumulative dist
std::vector<double> photovalue; //array to store values of cumulative dist
tk::spline photo_cross_table;
struct f_params {double a; tk::spline b; double c;};

void SetTable(tk::spline &f_cross,std::vector<double> &x, std::vector<double> &y,const string& filename) {

  ifstream file (filename.c_str());
  double d_x=0,d_y=0;
  string line;

  if (file.is_open()){
    while(getline (file,line) ){
      istringstream in(line);
      in>>d_x;
      in>>d_y;
      x.push_back(d_x);
      y.push_back(d_y);
    }
  }
  else cout << "File could not be opened" << endl;
  f_cross.set_points(x,y);

  return;
}

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
  photo_cross_table.set_points(photoenergy,photovalue);

  return;
}


double photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  if( x>=10.6406 && x<500 ) f = f_cross(x);
  else if( x >= 500 && x < 3206 ) f = 1.346*pow(500/x,2.54);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}


double LAP_photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  if( x>=15.75 && x<=245 ) f = f_cross(x);
  else if( x > 245 && x < 3206 ) f = 4.51*pow(248/x,2.29);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}


double photo_cross_interp(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  double emin = photoenergy.front();
  double emax = photoenergy.back();
  double f=0.;//function value
  if( x>=emin && x<=emax ) f = photo_cross_table(x);
  else f = 0;

  return f;

}



double dipole_oscill(double x, void *p) {
  double f = 0.;
  double z = 18; //atomic number
  double scale = 1e-18; //photo cross section is scaled down by 1e-18
  f = (LAP_photo_cross(x,p) * scale)/(1.097e-16 * z);

  return f;
}

double log_dipole_oscill(double x, void *p) {
  double f = 0.;
  double z = 18; //atomic number
  double scale = 1e-18; //photo cross section is scaled down by 1e-18
  f = log(x)*(LAP_photo_cross(x,p) * scale)/(1.097e-16 * z);

  return f;
}

int main(){

  ofstream output;
  std::vector<double> energy,cross; //vectors to store the table data
  tk::spline f_cross;

  //SetTable(f_cross,energy,cross,"./argon10eV_500eV.dat");
  SetTable(f_cross,energy,cross,"./argon_west_test.dat");
  //  SetPhotoCross("./argon_west.dat");
  //   SetPhotoCross("blum_rold_argon_photo.dat");
  
  int num_points = 5e2;

  double units_cross = 1e-18;  // [cm^2]
  double density = 1.662e-3;   // [g/cm^3]
  double molarmass = 39.948;   // [g/mol]
  //ASSUMPTION!!!!
  //we assume that the units of cross section, sigma, are given in cm^2
  //units of cross seciton read in are not SI units
  //we have to incorperate this scale somwhere
  //I chose to incorperate it into avogadros number
  //since avogadros number is so large to prevent precision error
  double avogadro = 6.022e23*units_cross; // [atoms/mol]*scalefactor
  double atom_cm3 = (density / molarmass) * avogadro; //[atoms/cm^3] or N in the equation

  double e = 0;//not used in this test
  struct f_params alpha = {e,f_cross,atom_cm3};
  
  double result, error;
  double expected = .001;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  
  F.function = &dipole_oscill;
  // F.function = &dipole_oscill;
  F.params = &alpha;
  gsl_integration_qag (&F, 0, 1e5, 0, 1e-3,1000,2,w, &result, &error); 
  
  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);
  
  gsl_integration_workspace_free (w);

  void *g = &alpha;
  int npoints = 1e5;
  double energystep = .1;
  TGraph interp = TGraph(npoints);

  for(int i=0;i<npoints;++i){
    double energy = energystep*i;
    double sigma = LAP_photo_cross(energy,g);
    //    if(energy>200 && energy<230)cout<<energy<<" "<<sigma<<endl;
    interp.SetPoint(i,energy,sigma);
  }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->SetLogx();
  c1->SetLogy();
  interp.Draw();
  c1->SaveAs("interp_cross.png");

  

  //fill the TGraphs
  void * p = &alpha;

  
  return 0;
  
}



