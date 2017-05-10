#include <gsl/gsl_integration.h>
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

double re_epsilon(double, void *);
double photo_cross(double, void *);
double photo_cross_interp(double, void *);
double LAP_photo_cross(double , void *);

double emin=0,emax=0;
std::vector<double> photoenergy; //array to store values of cumulative dist
std::vector<double> photovalue; //array to store values of cumulative dist
struct f_params {double a; gsl_spline  *b; double c;};

gsl_interp_accel *acc = gsl_interp_accel_alloc();
gsl_spline *photo_cross_table;

void SetTable(gsl_spline *f_cross,std::vector<double> &x, std::vector<double> &y,const string& filename) {

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
      cout<<d_energy<<" "<<d_value<<endl;
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  double *x_a = photoenergy.data();
  double *y_a = photovalue.data();

  int size_a = photoenergy.size();
  emin = photoenergy.front();
  emax = photoenergy.back();

  photo_cross_table = gsl_spline_alloc(gsl_interp_akima,size_a);
  gsl_spline_init(photo_cross_table,x_a,y_a,size_a);

  return;
}


double integrand(double x, void *p) {
  //def of QAWC integration in gsl library is
  //I = \int_a^b dx f(x) / (x - c)
  //we want to integrate I = \int_(0,infinity) dE' f(E') / (E'^2 - E^2)
  //which can also be simplified as f(E')/(E'-E)*(E'+E)
  //since we are integrating from 0 to infinity in E'
  //we can put f(E')/(E'+E)=f"(E')
  //thus the new integral is I = dE' f"(E')/(E'-E)
  //which is the same as the gsl format
  //with the singularity at (E'-E)

  struct f_params * params = (struct f_params *)p;
  double e = (params->a);//energy E in the integrand 

  double f =  (x * re_epsilon(x,p))/(x + e);
  f = f*(2/3.1415);

  return f;
}

double integrand3(double x, void *p) {

  struct f_params * params = (struct f_params *)p;
  double e = (params->a);//energy E in the integrand 
  
  double f = (e*re_epsilon(x,p)-e*re_epsilon(e,p));
  f = f/(pow(x,2)-pow(e,2));
  f = f*(2/3.1415);
  f = -f;
  
  return f;
}

double re_epsilon(double x, void * p) {
  double f=0.;//function value

  if( x>=emin && x<=emax ) f =  gsl_spline_eval(photo_cross_table,x,acc);
  else f = 0;

  return f;
}

int main(){

  ofstream output;
  std::vector<double> x,y; //vectors to store the table data
  gsl_spline *f_cross;

  //      SetPhotoCross("./bichsel_Re.dat");
      SetPhotoCross("./myCH4.dat");

  int num_points = 9e2;
  TGraph *imaginary = new TGraph(num_points);

  double units_cross = 1e-18;  // [cm^2]
  //  double density = 1.662e-3;   // [g/cm^3]
  //  double density = 1.783e-3;   // [g/cm^3]
  //  double molarmass = 39.948;   // [g/mol]
  double density = .6668e-3;   // [g/cm^3]
  double molarmass = 16.04246;   // [g/mol]
  //
  //ASSUMPTION!!!!
  //we assume that the units of cross section, sigma, are given in cm^2
  //units of cross seciton read in are not SI units
  //we have to incorperate this scale somwhere
  //I chose to incorperate it into avogadros number
  //since avogadros number is so large to prevent precision error
  double avogadro = 6.022e23*units_cross; // [atoms/mol]*scalefactor
  double atom_cm3 = (density / molarmass) * avogadro; //[atoms/cm^3] or N in the equation



  output.open("img_fromReal.dat");
  if (output.is_open()){
    output<<"Energy[eV]"<<"\t"<<"Img"<<endl;
   
    double energy_step = 1.e2/num_points;
    for(int i=1;i<num_points;++i){

      double e = i * energy_step; //[eV] value to evaluate epsilon(e) at E in the

      cout<<"Energy step is "<<e<<endl;
      struct f_params alpha = {e,f_cross,atom_cm3};
      
      double result, error;
      double expected = .001;

      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

      gsl_function F;
      //      F.function = &dipole_oscill;
      //      F.function = &integrand;
      F.function = &integrand3;
      F.params = &alpha;

      gsl_integration_qag (&F, .001, 1e4, 0, 1e-3,1000,6,w, &result, &error); 
      //      //      gsl_integration_qawc (&F, .1, 1e5, e, 0 , 1e-3, 1000,w, &result, &error); 
      /*
      printf ("result          = % .18f\n", result);
      printf ("exact result    = % .18f\n", expected);
      printf ("estimated error = % .18f\n", error);
      printf ("actual error    = % .18f\n", result - expected);
      printf ("intervals       = %zu\n", w->size);
      */
      gsl_integration_workspace_free (w);

      //fill the TGraphs

      double re_value = (result);
      imaginary->SetPoint(i,e,result);
      output<<e<<"\t"<<re_value<<endl;

    }
  }

  TCanvas *c2 = new TCanvas(1);
  c2->cd();
   imaginary->GetXaxis()->SetRangeUser(8,100);
   imaginary->GetYaxis()->SetRangeUser(-.0025,.0015);
  imaginary->SetTitle("#delta#[]{E} = Re[#epsilon#(){E}-1]");
  imaginary->GetYaxis()->SetTitle("#delta[E]*10^{4}");
  imaginary->GetYaxis()->CenterTitle();
  imaginary->GetXaxis()->SetTitle("E [eV]");
  imaginary->GetXaxis()->CenterTitle();
  imaginary->Draw();
  //  c2->SetLogy();
  c2->SetLogx();
  c2->SaveAs("imaginary_fromReal.png");

  return 0;

  }



