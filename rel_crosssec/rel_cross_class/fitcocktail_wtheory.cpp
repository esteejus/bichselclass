#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <functional>
#include <gsl/gsl_integration.h>
#include "TAxis.h"
#include "TString.h"
#include "TH1.h"
#include <algorithm>
#include <numeric>

using namespace std;

vector<double> data,expected;

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

double calcInt(double x1, double x2,vector<double> &a, vector<double> &b){

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
    cout<<x<<" "<<f<<endl;
    return f;};

  std::function<double(double)> F2(integral);
  gsl_function_pp F2_2(F2);
  gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_integration_qag (F_2, x1, x2, 0, 1e-3, 10000,6,w, &result, &error);    

  gsl_integration_workspace_free (w);
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);

  return result;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag){



  return;
}

vector<double> find_min(vector<double>dep, vector<double> ind, vector<string> dep_name, vector<string> ind_name, int npar){

  // Initialize minuit, set initial values etc. of parameters.
  //  const int npar = 2;              // the number of parameters
  vector <double> output ={};
  TMinuit minuit(npar);
  //  minuit.SetPrintLevel(-1);
  minuit.SetFCN(fcn);

  double par[npar];               // the start values
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];

  par[0] = 10;            // a guess at the true value.
  stepSize[0] = .1;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";
  
  par[1] = 22;            // a guess at the true value.
  stepSize[1] = .1;       // take e.g. 0.1 of start value
  minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[1] = 0;
  parName[1] = "scale";

  for (int i=0; i<npar; i++){
    minuit.DefineParameter(i, parName[i].c_str(), 
			   par[i], stepSize[i], minVal[i], maxVal[i]);
  }
  
  // Do the minimization!
  
  minuit.Migrad();       // Minuit's best minimization algorithm
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
    minuit.GetParameter(i,outpar[i],err[i]);
    output.push_back(outpar[i]);
  }

  cout << endl << endl << endl;
  cout << "   "<<parName[0]<<": " << outpar[0] << " +/- " << err[0] << endl;
  cout << "   "<<parName[1]<<": " << outpar[1] << " +/- " << err[1] << endl;


  return output;

}

void getExpected(vector<double> &a, vector<double> &b,TH1D *data, double par[]){

  double offset = par[0];
  double slope = par[1];
  double counts = par[2];

  expected_vec.clear();
  int npoints = data -> GetNbinsX();
  for(int i = 1;i <= npoints; i++)
    {
      //x1 & x2 are in units of data
      // data = theory * scale + offset
      //theory = (data - offset)/scale
      double x1 = data->GetBinLowEdge(i);
      double x2 = x1 + data -> GetBinWidth(i);
      x1 = (x1 - offset)/slope;
      x2 = (x2 - offset)/slope;

      double integral = calcInt(x1,x2,a,b);
      integral *= counts;
      expected_vec.push_back(integral);
    }

  return;
}

TGraph * GraphFunc(vector<double> &a, vector<double> &b){
  int size_x = a.size();
  int npoints = 2*a.size();
  TGraph * f = new TGraph(npoints);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
  
  double x1 = a.front(), x2 = a.back();
  if(x1 < 1e-3) x1 = .1;//using log cant have nan errors
  double log_energy_step = (log(x2)-log(x1))/npoints;
  for(int i = 0 ;i < npoints; i++)
    {
      double delta = i*log_energy_step + log(x1) ;
      delta = exp(delta);
      double value = 0;
      if(delta < a.back() && delta >a.front())
	value = gsl_spline_eval(dist_table,delta,acc);

      //      if(value<0)cout<<"Value interpolat is neg"<<delta<<" "<<value<<endl;
      f->SetPoint(i,delta,value);
    }
  
  return f;
}

int main(){

  TFile *f = new TFile("1700MeV_cocktail.root","OPEN");
  TH1D *data = (TH1D*)f ->Get("dedxdist");

  TFile *g = new TFile("cocktail300.root","OPEN");
  TH1D *p_dist = (TH1D*)g ->Get("p_c_dist");
  TH1D *d_dist = (TH1D*)g ->Get("d_c_dist");
  TH1D *t_dist = (TH1D*)g ->Get("t_c_dist");

  //  double slope  = 19.2335;
  double slope  = 26;
  double offset = -4;
  //data = theory*slope + offset
  // theory = (data - offset)/slope

    int binl = data->GetXaxis()->FindBin(20);
    int binh = data->GetXaxis()->FindBin(50);
    double entries = data->Integral(binl,binh,"width");
    cout<<entries<<endl;

  vector<double> p_x, p_y, d_x,d_y, t_x,t_y;
  for(int i = 1 ;i < p_dist->GetNbinsX(); i++)
    {
      double x = p_dist->GetBinCenter(i);
      double y = p_dist->GetBinContent(i);
      x /= 1e3;
      x = x*slope + offset;
      if(y>1e-5)
	{
	  p_x.push_back(x);
	  p_y.push_back(y);
	}

    }

  double scale = calcInt(p_x,p_y);
  cout<<"new scale"<<scale<<endl;
  for(int i = 0 ;i < p_y.size(); i++)
    p_y.at(i) = p_y.at(i)*(entries/scale);

  /*
  for(int i = 1 ;i < t_dist->GetNbinsX(); i++)
    {
      double x = t_dist->GetBinCenter(i);
      double y = t_dist->GetBinContent(i);
      x /= 1e3;
      x = x*slope + offset;
      if(y<1e-4)
	{
	  t_x.push_back(x);
	  t_y.push_back(y);
	}
    }

    for(int i = 1 ;i < d_dist->GetNbinsX(); i++)
    {
      double x = d_dist->GetBinCenter(i);
      double y = d_dist->GetBinContent(i);
      x /= 1e3;
      x = x*slope + offset;
      if(y<1e-4)
	{
	  d_x.push_back(x);
	  d_y.push_back(y);
	}
    }
  */
  //  double scale = calcInt(d_x,d_y);
    //  cout<<"new scale"<<scale<<endl;
  //  for(int i = 0 ;i < d_y.size(); i++)
  //    d_y.at(i) = d_y.at(i)*(entries/scale);

  //TGraph *p_scaled = GraphFunc(p_x,p_y);
    TH1D *p_hist = getHist(p_x,p_y,16,offset,slope,1);
    //  TH1D *d_hist = getHist(d_x,d_y,13,offset,slope,1);

  //    TGraph *t_scaled = GraphFunc(t_x,t_y);

    TCanvas *c1 = new TCanvas("","",1);
    c1->SetLogx();
        data -> GetXaxis() ->SetRangeUser(20,50);
        data -> Draw();

    //    d_scaled -> SetLineColor(2);
    //    t_scaled -> SetLineColor(4);
    //    p_scaled->GetXaxis()->SetRangeUser(0.
        p_hist->SetLineColor(2);
        p_hist -> Draw("same");
    //    d_hist -> Draw("same");
    //    d_scaled -> Draw("same");
    //    t_scaled -> Draw("same");

    c1->SaveAs("fitcocktail.png");

  return 0;
}
