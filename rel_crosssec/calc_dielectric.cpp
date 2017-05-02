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



/*
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
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  double z = 18;
  //    double f = ( (2/3.1415) * x * im_epsilon(x,p))/(x + e);
  //    double f = ((n*2*hbar_c)/3.1415)*photo_cross(x,p)/(x+e);
    double f = ((n*2*hbar_c)/(3.1415))*photo_cross_interp(x,p)/(x+e);
  //  double f = ((n*2*hbar_c)/3.1415)*photo_cross(sqrt(x),p)/(2*sqrt(x));
  return f;
}

*/
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
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]

  //  double f = ((n*2*hbar_c)/(3.1415*z))*photo_cross_interp(x,p)/(x+e);
  double f =  (x * im_epsilon(x,p))/(x + e);
  f = f*(2/3.1415);

  return f;
}


double integrand2(double x, void *p) {
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
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  double z = 18;
  double f = ((n*2*hbar_c*x)/(3.1415*z))*photo_cross(x,p)/(pow(x,2)-pow(e,2));

  return f;
}


/*
double integrand3(double x, void *p) {
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
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  double z = 18;
  //  double f = ((n*2*hbar_c)/(3.1415*z))*(x*photo_cross(x,p)-e*photo_cross(e,p))/(pow(x,2)-pow(e,2));
  
  //    double f = (x*photo_cross_interp(x,p)-e*photo_cross_interp(e,p));
    double f = (x*im_epsilon(x,p)-e*im_epsilon(e,p));
  //  double f = (x*photo_cross(x,p)-e*photo_cross(e,p));
  //  double f = (x*LAP_photo_cross(x,p)-e*LAP_photo_cross(e,p));
  f= f/(pow(x,2)-pow(e,2));
  f=f*(2/3.1415);
    
  return f;
}
*/

double integrand3(double x, void *p) {
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
  double n = (params->c);// [atoms/cm^3]*scale
  double hbar_c = 1.97327e-5; // [eV *cm]
  
double f = (x*im_epsilon(x,p)-e*im_epsilon(e,p));
 f = f/(pow(x,2)-pow(e,2));
 f = f*(2/3.1415);
    
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
  //  f = coeff*LAP_photo_cross(x,p)/x;
  f = coeff*photo_cross_interp(x,p)/x;
  return f;
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
  if( x>=15.75 && x<=248 ) f = f_cross(x);
  else if( x > 248 && x < 3206 ) f = 4.51*pow(248/x,2.29);
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
  f = (photo_cross(x,p) * scale)/(1.097e-16 * z);

  return f;
}

double log_dipole_oscill(double x, void *p) {
  double f = 0.;
  double z = 18; //atomic number
  double scale = 1e-18; //photo cross section is scaled down by 1e-18
  f = log(x)*(photo_cross(x,p) * scale)/(1.097e-16 * z);

  return f;
}

int main(){

  ofstream output;
  std::vector<double> energy,cross; //vectors to store the table data
  tk::spline f_cross;

  //  SetTable(f_cross,energy,cross,"./argon10eV_500eV.dat");
  //  SetTable(f_cross,energy,cross,"./argon_west.dat");
  SetPhotoCross("./argon_west.dat");
  //  SetPhotoCross("blum_rold_argon_photo.dat");
  
  int num_points = 5e2;
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
      
      double e = i * energy_step; //[eV] value to evaluate epsilon(e) at E in the
      if(e<10)continue;
      cout<<"Energy step is "<<e<<endl;
      struct f_params alpha = {e,f_cross,atom_cm3};
      
      double result, error;
           double expected = .001;

      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

      gsl_function F;
            F.function = &integrand;
      //      F.function = &integrand3;
      F.params = &alpha;
      //gsl_integration_qag (&F, 0, 1e5, 0, 1e-3,1000,2,w, &result, &error); 
      gsl_integration_qawc (&F, 1, 1e5, e, 0 , 1e-3, 1000,w, &result, &error); 
      cout<<i<<endl;
      printf ("result          = % .18f\n", result);
      printf ("exact result    = % .18f\n", expected);
      printf ("estimated error = % .18f\n", error);
      printf ("actual error    = % .18f\n", result - expected);
      printf ("intervals       = %zu\n", w->size);

      gsl_integration_workspace_free (w);

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
    cout<<im_energystep<<endl;
      for(int i=1;i<=im_steps;++i){
	
	double im_energy = i * im_energystep;//point E to evaluate
	struct f_params alpha = {im_energy,f_cross,atom_cm3};
	void * p = &alpha;
	
	double im_value = im_epsilon(im_energy,p);
		cout<<im_value<<endl;
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



