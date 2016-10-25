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


double photo_cross(double, void *);
double dipole_oscill(double, void *);
struct f_params { tk::spline b; double c; double v;};

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

double cross_section(double x, void * p) {

  double f = 0.;
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b); //getting the cross

  double z = 19;       //atomic number
  double a = 40;       //atomic mass g/mol
  double mass = .511;  //mass of electron MeV/c^2
  double beta = (params->v);

  //for heavy particles Emax = 2mc^2beta^2 gamma^2
  //Emax = T/2 for electrons
  double e_max = 2 * 5.11e5 * ( pow(beta,2) /( 1-pow(beta,2) ) );
  //cross section is zero when E>Emax
  cout<<"emax is "<<e_max<<endl;

  if(x<e_max){
    cout<<"here?"<<endl;
    double coeff = (1.5354e5 * z / ( a * pow(beta,2) ));
    double ruth_cross = coeff * ((1 - ( pow(beta,2) * x / e_max))/pow(x,2)); 

    double result, error;
    double expected = 158.176;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &dipole_oscill;
    F.params = &params;
    cout<<"what is x "<<x<<endl;
    gsl_integration_qags (&F, 0, x, 0, 1e-6, 1000,w, &result, &error); 

    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals       = %zu\n", w->size);

    f = (x * dipole_oscill(x,p) * log( (2 * mass * pow(beta,2)) / x) ) + result;
    f *= ruth_cross;

    gsl_integration_workspace_free (w);
  }

  else f = 0;

  return f;

}

double photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  if( x>=10.6406 && x<=500 ) f = f_cross(x);
  else if( x > 500 && x < 3206 ) f = 1.346*pow(500/x,2.54);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}

double dipole_oscill(double x, void *p) {
  double f = 0.;
  double z = 18; //atomic number
  double scale = 1e-18; //photo cross section is scaled down by 1e-18
  f = log(x)*(photo_cross(x,p) * scale)/(1.097e-16 * z);
  //  f = (photo_cross(x,p) * scale)/(1.097e-16 * z);

  return f;
}

int main(){

  ofstream output;
  std::vector<double> energy,cross; //vectors to store the table data
  tk::spline f_cross;

  SetTable(f_cross,energy,cross,"./argon10eV_500eV.dat");

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
  double beta_gamma = .2;
  double beta = sqrt(pow(beta_gamma,2)/(1 + pow(beta_gamma,2)));

  struct f_params alpha = {f_cross,atom_cm3,beta};
  
    /*
    double result, error;
    double expected = 158.176;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &im_epsilon;
    F.params = &alpha;
    //  gsl_integration_qags (&F, 500, 10000, 0, 1e-6, 1000,w, &result, &error); 
    gsl_integration_qawc (&F, 1, 1e6, e, 0, 1e-3, 1000,w, &result, &error); 

    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

  }
}

    */

  int steps = 1e4;
  double energystep = 1.e4/steps;
  TGraph *crosssection = new TGraph(steps);
  void * p = &alpha;
  double f = 0.;
  
  for(int i=1;i<=steps;++i){
    double e_value = i * energystep;
    //    double value = cross_section(e_value,p);
    //    crosssection->SetPoint(i,e_value,value);

    f = 0.;

  double z = 19;       //atomic number
  double a = 40;       //atomic mass g/mol
  double mass = .511;  //mass of electron MeV/c^2

  //for heavy particles Emax = 2mc^2beta^2 gamma^2
  //Emax = T/2 for electrons
  double e_max = 2 * 5.11e5 * ( pow(beta,2) /( 1-pow(beta,2) ) );
  //cross section is zero when E>Emax

  double result, error;
  double expected = 0;

  if(e_value<e_max){

    double coeff = (1.5354e5 * z / ( a * pow(beta,2) ));
    double ruth_cross = coeff * ((1 - ( pow(beta,2) * e_value / e_max))/pow(e_value,2)); 

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &dipole_oscill;
    F.params = &alpha;
    cout<<"what is x "<<e_value<<endl;
    gsl_integration_qags (&F, 0, e_value, 0, 1e-3, 1000,w, &result, &error); 

    printf ("result          = % .18f\n", result);
    printf ("exact result    = % .18f\n", expected);
    printf ("estimated error = % .18f\n", error);
    printf ("actual error    = % .18f\n", result - expected);
    printf ("intervals       = %zu\n", w->size);

    f = (e_value * dipole_oscill(e_value,p) * log( (2 * mass * pow(beta,2)) / e_value) ) + result;
    f *= ruth_cross;
    gsl_integration_workspace_free (w);
  }

  else f = 0;

  crosssection->SetPoint(i,e_value,f);


  }
 
  
TCanvas *c1 = new TCanvas(1);
//crosssection->SetTitle("#delta#[]{E} = Im[#epsilon#(){E}]");
//crosssection->GetYaxis()->SetTitle("#delta#[]{E}*10^{4}");
crosssection->GetYaxis()->CenterTitle();
crosssection->GetXaxis()->SetTitle("E [eV]");
crosssection->GetXaxis()->CenterTitle();
 crosssection->GetXaxis()->SetRangeUser(10,500);
crosssection->Draw();
 
c1->SetLogx();
c1->SetLogy();
c1->SaveAs("crosssection.jpg");
  
return 0;

}



