#include <gsl/gsl_integration.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "spline.h"
#include <cmath>


using namespace std;

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

double integrand(double x, void * params) {
  double f = 0.;
  double e = *(double *) params;
  //def of QAWC integration in gsl library is
  //I = \int_a^b dx f(x) / (x - c)
  //we want to integrate I = \int_a^b dx f(E) / (E'^2 - E^2)
  //there fore x in this integrand function represents E'^2
  //we must therefore take sqrt of x to get E'
  //this is the new variable to fit the gsl format
  f = photo_cross(sqrt(x),params);

  return f;
}

double photo_cross(double x, void * params) {

  tk::spline f_cross = *(tk::spline *) params;
  double f=0.;//function value
  if( x>=10.6406 && x<=500 ) f = f_cross(x);
  else if( x > 500 && x < 3206 ) f = 1.346*pow(500/x,2.54);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}


/*
double photo_cross(double x, void *params){
  double f =0.;
  f = 1.346*pow(500/x,2.54);
  return f;
}
*/

int main(){

  std::vector<double> energy,cross; //vectors to store the table data
  tk::spline f_cross;
  SetTable(f_cross,energy,cross,"./argon10eV_500eV.dat");

gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, error;
  double expected = 158.176;

  gsl_function F;
  F.function = &photo_cross;
  F.params = &f_cross;

  gsl_integration_qags (&F, 500, 10000, 0, 1e-6, 1000,
                        w, &result, &error); 

  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free (w);



  
  return 0;
}



