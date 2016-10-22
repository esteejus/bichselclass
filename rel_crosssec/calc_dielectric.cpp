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

double integrand(double x, void *p) {
  struct f_params * params = (struct f_params *)p;
  double e = (params->a);//energy E in the integrand 
  double f = (2/3.1415) * x * im_epsilon(x,p)/(x + e);
    
  return f;
}
  
double im_epsilon(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  //  double e = (params->a);//energy E in the integrand 
  //n is usually scaled to account for the real values of cross section
  double n = (params->c);// [atoms/cm^3]*scale
  //  double b = (2/3.1415)*1.97327e-5; //[cm *eV] b as defined in Lapique B pg.18 
  double hbar_c = 1.97327e-5; // [eV *cm]
  double coeff = hbar_c*n;

  double f = 0.;
  //def of QAWC integration in gsl library is
  //I = \int_a^b dx f(x) / (x - c)
  //we want to integrate I = \int_(0,infinity) dE' f(E') / (E'^2 - E^2)
  //which can also be simplified as f(E')/(E'-E)*(E'+E)
  //since we are integrating from 0 to infinity in E'
  //we can put f(E')/(E'+E)=f"(E')
  //thus the new integral is I = dE' f"(E')/(E'-E)
  //which is the same as the gsl format
  //with the singularity at (E'-E)

  f = coeff*photo_cross(x,p)/x;

  return f;
}

double photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);

  double f=0.;//function value
  if( x>=10.6406 && x<=500 ) f = f_cross(x);
  else if( x > 500 && x < 3206 ) f = 1.346*pow(500/x,2.54);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}

int main(){

  std::vector<double> energy,cross; //vectors to store the table data
  tk::spline f_cross;
  SetTable(f_cross,energy,cross,"./argon10eV_500eV.dat");

  int num_points = 1e3;
  TGraph *imaginary = new TGraph(num_points);
  TGraph *real = new TGraph(num_points);//actually Re[e]-1


  for(int i=5;i<num_points;++i){
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
    double result, error;
    double expected = 158.176;
    
    double units_cross = 1e-18;  // [cm^2]
    double density = 1.784e-3;   // [g/cm^3]
    double molarmass = 39.948;   // [g/mol]
    //ASSUMPTION!!!!
    //we assume that the units of cross section, sigma, are given in cm^2
    //units of cross seciton read in are not SI units
    //we have to incorperate this scale somwhere
    //I chose to incorperate it into avogadros number
    //since avogadros number is so large to prevent precision error
    double avogadro = 6.0022e23*units_cross; // [atoms/mol]*scalefactor
    double e = i; //[eV] value to evaluate epsilon(e) at E in the integrand
    double atom_cm3 = (density / molarmass) * avogadro; //[atoms/cm^3] or N in the equation

    cout <<" Number of atoms/cm^3 is "<<atom_cm3<<endl;
    struct f_params alpha = {e,f_cross,atom_cm3};

    gsl_function F;
    F.function = &im_epsilon;
    F.params = &alpha;

    //  gsl_integration_qags (&F, 500, 10000, 0, 1e-6, 1000,w, &result, &error); 
    gsl_integration_qawc (&F, 4, 5000, e, 0, 1e-3, 1000,w, &result, &error); 

    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

    //fill the TGraphs
    void * p = &alpha;
    imaginary->SetPoint(i,i,im_epsilon(i,p)/1e-4);
    real->SetPoint(i,e,result/1e-4);
  }


  TCanvas *c1 = new TCanvas(1);
 imaginary->SetTitle("#delta#[]{E} = Im[#epsilon#(){E}]");
 imaginary->GetYaxis()->SetTitle("#delta#[]{E}*10^{4}");
 imaginary->GetYaxis()->CenterTitle();
 imaginary->GetXaxis()->SetTitle("E [eV]");
 imaginary->GetXaxis()->CenterTitle();
 imaginary->Draw();
 
 c1->SetLogx();
 c1->SaveAs("imaginary.jpg");
 
 TCanvas *c2 = new TCanvas(1);
 c2->cd();
 //  real->GetXaxis()->SetRangeUser(100,2000);
 // real->GetYaxis()->SetRangeUser(-.1,.06);
 real->SetTitle("#delta#[]{E} = Re[#epsilon#(){E}-1]");
 real->GetYaxis()->SetTitle("#delta[E]*10^{4}");
 real->GetYaxis()->CenterTitle();
 real->GetXaxis()->SetTitle("E [eV]");
 real->GetXaxis()->CenterTitle();
 real->Draw();
 c2->SetLogx();
  c2->SaveAs("real.jpg");
  
  return 0;
}



