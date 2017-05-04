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
double LAP_photo_cross(double , void *);
double discrete_sum(double);
  
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


double integrand(double x, void *p) {

  struct f_params * params = (struct f_params *)p;
  double e = (params->a);//energy E in the integrand 
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  
  double f = (x*im_epsilon(x,p)-e*im_epsilon(e,p));
  f = f/(pow(x,2)-pow(e,2));
  f = f*(2/3.1415);

  //  f = 0;
  //  double discrete = discrete_sum(e);    
  //  cout<<"Discrete value is "<<discrete<<endl;
  //  discrete = discrete*n;
  //  f += discrete;
  
  return f;
}

double integrand_c(double x, void *p) {

  struct f_params * params = (struct f_params *)p;
  double e = (params->a);//energy E in the integrand 
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  
  double f = x*im_epsilon(x,p);
  f = f/(x + e);
  f = f*(2/3.1415);

  //  double discrete = discrete_sum(e);    
  //  cout<<"Discrete value is "<<discrete<<endl;
  //  discrete = discrete*n;
  //  f += discrete;
  
  return f;
}

double im_epsilon(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  //n is usually scaled to account for the real values of cross section
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  double z      = 1;
  double coeff = (hbar_c * n)/z;
  double f = 0.;
  f = coeff*LAP_photo_cross(x,p)/x;

  return f;
}

double LAP_photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  if( x>=15.75 && x<245 ) f = f_cross(x);
  else if( x >= 245 && x < 3206 ) f = 4.51*pow(248/x,2.29);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;


  return f;

}

double discrete_sum(double x){
  double f =0;
  double z = 18;
  vector<double> f_strength = {59     ,228   ,28    ,13     ,93     , 107   };//not scaled/normalized by z
  vector<double> energy     = {11.6214,11.828,14.122,14.2551,14.1525,14.3035};

  double hbar_c = 1.97327e-5; // [eV *cm]
  double coeff = (2/3.1415926)*hbar_c*1.097e-16;//A in LAP paper

  for(int i = 0; i < f_strength.size(); ++i) {
    f += (f_strength.at(i)/z)/(pow(x,2)-pow(energy.at(i),2));
  }

  return (f*coeff*1e-3);  //paper quotes X 10e3
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

  //    SetTable(f_cross,energy,cross,"./blum_rold_argon_photo.dat");
  //        SetTable(f_cross,energy,cross,"./argon10eV_500eV.dat");
  SetTable(f_cross,energy,cross,"./argon_west_test.dat");
  //  SetPhotoCross("./argon_west.dat");
  //  SetPhotoCross("blum_rold_argon_photo.dat");
  
  int num_points = 1e4;
  TGraph *imaginary = new TGraph(num_points);
  TGraph *real = new TGraph(num_points);//actually Re[e]-1

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



  output.open("real.dat");
  if (output.is_open()){
    output<<"Data table for Real componets Re[epsilon] for Argon gas"<<endl;
    output<<"Re values are multiplied by 1e4 to better display the scale"<<endl;
    output<<"Energy[eV]"<<"\t"<<"Re*1e4"<<endl;
   
    double energy_step = 2.e2/num_points;
    for(int i=1;i<num_points;++i){
      
      double e = i * energy_step; //[eV] value to evaluate epsilon(e) at E in
      if(e<15.9)continue;
      cout<<"Energy step is "<<e<<endl;
      struct f_params alpha = {e,f_cross,atom_cm3};
      
      double result, error;
           double expected = .001;

      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

      gsl_function F;
      F.function = &integrand_c;
      F.params = &alpha;
      gsl_integration_qawc (&F, 15.7, 1e5, e, 0 , 1e-4, 1000,w, &result, &error); 
      //      gsl_integration_qag (&F, 15.7, 1e5, 0, 1e-3,1000,2,w, &result, &error); 
      printf ("result          = % .18f\n", result);
      printf ("exact result    = % .18f\n", expected);
      printf ("estimated error = % .18f\n", error);
      printf ("actual error    = % .18f\n", result - expected);
      printf ("intervals       = %zu\n", w->size);

      gsl_integration_workspace_free (w);

      double discrete = discrete_sum(e)*atom_cm3;
      result += discrete;
      cout<<"discrete value is "<<discrete<<endl;
      //fill the TGraphs
      void * p = &alpha;

      double re_value = (result);
      //fill TGraph
      //      real->SetPoint(i,e,result/1e-4);
      real->SetPoint(i,e,result);
      //output to .dat file
      output<<e<<"\t"<<re_value<<endl;


    }
  }

  ofstream im_output;
  im_output.open("imaginary.dat");

  if (im_output.is_open()){
    int im_steps = 1e6;
    double im_energystep = 10000./im_steps;
    //    cout<<im_energystep<<endl;
      for(int i=1;i<=im_steps;++i){
	
	double im_energy = i * im_energystep;//point E to evaluate
	struct f_params alpha = {im_energy,f_cross,atom_cm3};
	void * p = &alpha;
	
	double im_value = im_epsilon(im_energy,p);
	//		cout<<im_value<<endl;
	imaginary->SetPoint(i,im_energy,im_value);
	im_output<<im_energy<<"\t"<<im_value<<endl;
      }
  }

  cout<<"NUM of atoms/cm^3 is "<<atom_cm3<<endl;
  
  TCanvas *c1 = new TCanvas(1);
  imaginary->SetTitle("#delta#[]{E} = Im[#epsilon#(){E}]");
  imaginary->GetYaxis()->SetTitle("#delta#[]{E}*10^{4}");
  imaginary->GetYaxis()->CenterTitle();
  imaginary->GetXaxis()->SetTitle("E [eV]");
  imaginary->GetXaxis()->CenterTitle();
  //  imaginary->GetXaxis()->SetRangeUser(10,20);
    imaginary->GetXaxis()->SetRangeUser(10,1e3);
  imaginary->Draw();
 
  c1->SetLogx();
  c1->SetLogy();
  c1->SaveAs("imaginary.png");
 
  TCanvas *c2 = new TCanvas(1);
  c2->cd();
   real->GetXaxis()->SetRangeUser(.1,10000);
  // real->GetYaxis()->SetRangeUser(-.1,.06);
  real->SetTitle("#delta#[]{E} = Re[#epsilon#(){E}-1]");
  real->GetYaxis()->SetTitle("#delta[E]*10^{4}");
  real->GetYaxis()->CenterTitle();
  real->GetXaxis()->SetTitle("E [eV]");
  real->GetXaxis()->CenterTitle();
  real->Draw();
  c2->SetLogx();
  c2->SaveAs("real.png");
  
  return 0;

  }



