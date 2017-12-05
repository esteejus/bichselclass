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


TGraph * GraphFunc (vector<double> &a, vector<double> &b){
  int npoints = 1e2;
  int size_x = a.size();
  TGraph * f = new TGraph(npoints);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);

  
  double x1 = a.front(), x2 = a.back();
  double stepsize = (x2-x1)/npoints;

  for(int i = 0 ;i < npoints; i++)
    {
      double delta = stepsize*i + x1;
      double value = 0;
      if(delta < a.back() && delta >a.front())
	value = gsl_spline_eval(dist_table,delta,acc);

      //      if(value<0)cout<<"Value interpolat is neg"<<delta<<" "<<value<<endl;
      f->SetPoint(i,delta,value);
    }
  
  return f;
}

double calcInt(vector<double> &a, vector<double> &b){
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
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_integration_qag (F_2, 0, 1e3, 0, 1e-2, 1000,1,w, &result, &error);    

  gsl_integration_workspace_free (w);
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);

  return result;
}

void ConvSelf(vector<double> &a, vector<double> &b){
  int size_x = a.size();
      gsl_interp_accel *acc = gsl_interp_accel_alloc();
      gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
      gsl_spline_init(dist_table,a.data(),b.data(),size_x);

  for(int iEnergy = 0 ;iEnergy < size_x; iEnergy++){
    double result = 0, error = 0;
    double delta = a.at(iEnergy) ;

    auto integral = [&](double x)->double{
      double f = 0, f_d = 0;
      if(x >= a.front() && x <= a.back())
	f = gsl_spline_eval(dist_table,x,acc);
      else
	f = 0;
      if((delta-x) >= a.front() && (delta-x) <= a.back())
	f_d = gsl_spline_eval(dist_table,delta-x,acc);
      else
	f_d = 0;
      if(f_d*f<0)cout<<"f "<<"f_d "<<f<<" "<<f_d<<" "<<x<<" "<<delta-x<<endl;
      return f*f_d;};

    std::function<double(double)> F2(integral);
    gsl_function_pp F2_2(F2);
    gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag (F_2, 0, delta, 0, 1e-3, 1000,1,w, &result, &error);    
    b.at(iEnergy)=result;
    //    cout<<"result"<<result<<endl;
    gsl_integration_workspace_free (w);
  }
    gsl_interp_accel_free(acc);
    gsl_spline_free(dist_table);
      
}


int main(){

  TF1 * gaus = new TF1("gaus","[0]*TMath::Gaus(x,[1],[2])",-20,20);
  gaus -> SetParameters(1,10,2);
  double integral = gaus-> Integral(-20,20);
  gaus ->SetParameter(0,1./integral);
  integral = gaus-> Integral(-20,20);
  cout<<"integral "<<integral<<endl;
 
  vector<double> x,y;

  double x1 = -20, x2 = 80;
  int npoints = 1e3;
  double stepsize = (x2-x1)/npoints;

  for(int i = 0 ;i < npoints; i++)
    {
      double delta = stepsize*i + x1;
      double value = gaus->Eval(delta);
      x.push_back(delta);
      y.push_back(value);
    }

  TGraph *f1 = GraphFunc(x,y);
  cout<<"F1 int "<<calcInt(x,y)<<endl;
  ConvSelf(x,y);
  TGraph *f2 = GraphFunc(x,y);   
  cout<<"F2 int "<<calcInt(x,y)<<endl;
  ConvSelf(x,y);
  TGraph *f3 = GraphFunc(x,y);   
  cout<<"F3 int "<<calcInt(x,y)<<endl;
 
  TCanvas *c1 = new TCanvas("","",1);
  c1->SetLogy();
  f1->Draw();
  f2->SetLineColor(2);
  f2->Draw("same");
  f3->SetLineColor(4);
  f3->Draw("same");

  c1 -> SaveAs("convfunc.png");

  return 0;

}
   
