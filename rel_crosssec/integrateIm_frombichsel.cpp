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

double im_epsilon(double, void *);
double photo_cross(double, void *);
double photo_cross_interp(double, void *);
double LAP_photo_cross(double , void *);

double emin=0,emax=0;
std::vector<double> photoenergy; //array to store values of cumulative dist
std::vector<double> photovalue; //array to store values of cumulative dist

gsl_interp_accel *acc = gsl_interp_accel_alloc();
gsl_spline *photo_cross_table;

struct f_params {double a; gsl_spline *b; double c;};

void SetTable(gsl_spline *f_cross, vector<double> &x , vector<double> &y,const string& filename) {

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

  /*
double *x_a = &x[0];
  double *y_a = &y[0];

  int size_a = x.size();
  f_cross = gsl_spline_alloc(gsl_interp_cspline,size_a);
  gsl_spline_init(f_cross,x_a,y_a,size_a);
  */
  
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

  double *x_a = photoenergy.data();
  double *y_a = photovalue.data();

  int size_a = photoenergy.size();
  emin = photoenergy.front();
  emax = photoenergy.back();

  photo_cross_table = gsl_spline_alloc(gsl_interp_akima,size_a);
  gsl_spline_init(photo_cross_table,x_a,y_a,size_a);
  
  return;
}


double photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  gsl_spline *f_cross = (params->b);
  //  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  /*
  if( x>=10.6406 && x<500 ) f = f_cross(x);
  else if( x >= 500 && x < 3206 ) f = 1.346*pow(500/x,2.54);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;
  */

  if( x>=10.64058 && x<500 ) f = gsl_spline_eval(f_cross,x,acc);
  else if( x >= 500 && x < 3206 ) f = 4.51*pow(248/x,2.29);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

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
  f = photo_cross_interp(x,p);

  return f;
}

double LAP_photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  gsl_spline *f_cross = (params->b);
  //  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  if( x>=15.75 && x<248 ) f = gsl_spline_eval(f_cross,x,acc);
  else if( x >= 248 && x < 3206 ) f = 4.51*pow(248/x,2.29);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}


double photo_cross_interp(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  double f=0.;//function value

  if( x>=emin && x<=emax ) f = gsl_spline_eval(photo_cross_table,x,acc);
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
  std::vector<double> x,y; //vectors to store the table data

  gsl_spline *f_cross;
  SetPhotoCross("./bichsel_Re.dat");


  double *x_a = x.data();
  double *y_a = y.data();
  int size_a = x.size();
  //  f_cross = gsl_spline_alloc(gsl_interp_cspline,size_a);
  //  f_cross = gsl_spline_alloc(gsl_interp_akima,size_a);
  //  gsl_spline_init(f_cross,x_a,y_a,size_a);

  double units_cross = 1e-18;  // [cm^2]
  double density = .668e-3;   // [g/cm^3]
  double molarmass = 16.04246;   // [g/mol]
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
  void *g = &alpha;

  int num_int_points = 1e4;
  int j_integral = 1e5;

  double e_max = 9.9e4;
  double e_min = 8.;
  double h = (e_max - e_min)/j_integral;
  
  TGraph real = TGraph(num_int_points);

  output.open("real.dat");

  for(int i = 0; i < num_int_points; ++i)
    {
      if(i%1000==0)cout<<i<<endl;
      double sum = 0;
      double e_i = e_min + i*h;
      //      cout<<"Energy at "<<e_i<<endl;
      for(int j = 0;j < j_integral;++j)
	{
	  double e_j = e_min + j*h;
	  double e_j1 = e_min + (j+1)*h;
	  double im = im_epsilon(e_j,g);
	  double f = im/(e_j-e_i);
	  f += im/(e_j+e_i);
	  f *= .5;
	  f = -f;
	  //	  cout<<"img is "<<im<<" "<<e_j<<endl;
	  //cout<<"f is "<<f<<" "<<e_i<<" "<<e_j<<endl;
	  if(i%2 == 0)
	    {
	      if(!(j%2 == 0))
		{
		  sum += f;
		  // cout<<"here i should be even and j odd"<<i<< " j "<<j<<endl;
		}
	    }
	  else
	    {
	      if(j%2 == 0)
		{
		sum += f;
		//		cout<<"here i should be odd  and j even"<<i<< " j "<<j<<endl;
		}
	    }
	}

      sum *= (2/3.1415296)*(2*h);
      cout<<"Sum is "<<e_i<<" "<<sum<<endl;
      real.SetPoint(i,e_i,sum);
      output<<e_i<<"\t"<<sum<<endl;
    }


  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->SetLogx();
  c1->SetLogy();
  real.GetXaxis()->SetRangeUser(8,100);
  real.Draw();
  c1->SaveAs("real_maclaurins_direct.png");

  
  
  
  return 0;
  
}



