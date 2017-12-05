#include <gsl/gsl_spline.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <functional>
#include <gsl/gsl_integration.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include <algorithm>

using namespace std;

class gsl_function_pp : public gsl_function
 {
    private:
   std::function<double(double)> _func;
    static double invoke(double x, void *params) {
     return static_cast<gsl_function_pp*>(params)->_func(x);
   }

    public:
 gsl_function_pp(std::function<double(double)> const& func) : _func(func){
      function=&gsl_function_pp::invoke;
      params=this;
    }    
};

vector<double>photoenergy,photovalue;

void SetDataTable (const std::string& filename) {
  ifstream file (filename.c_str());
  double prev= 0;
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      if(d_energy<=prev)
      	cout<<" HEre "<<d_energy<<" "<<prev<<endl;

      //      cout<<d_energy<<" "<<d_value<<endl;
      photoenergy.push_back(d_energy/1000);
      photovalue.push_back(d_value);
      prev = d_energy;
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  return;
}

double calcInt(vector<double> &a, vector<double> &b){
  double x1 = a.front();
  double x2 = a.back();  
  int size_x = a.size();
  double result = 0, error = 0;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
      
  auto integral = [&](double x)->double{
    double f = 0;
    if(x >= a.front() && x <= a.back())
      f = gsl_spline_eval(dist_table,x,acc);
    else
      f = 0;

    return f;};

  std::function<double(double)> F2(integral);
  gsl_function_pp F2_2(F2);
  gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_integration_qag (F_2, x1, x2, 0, 1e-2, 10000,3,w, &result, &error);    

  gsl_integration_workspace_free (w);
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);

  return result;
}


int main(){

  SetDataTable("bichselStragCompare.dat");
  cout<<"Integral is "<<calcInt(photoenergy,photovalue);

  return 0;

}
