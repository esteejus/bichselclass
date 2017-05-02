#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include "spline.h"

#include <TMinuit.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include "TH2.h"
#include "TFile.h"
#include <TAxis.h>
#include <TLine.h>

using namespace std;


gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
gsl_spline *spline;// = gsl_spline_alloc (gsl_interp_cspline, size);

struct f_params {double m; double b;};//mx+b
tk::spline dedx_dist;
vector<double>t_dedx,t_value;//dedx value in MeV??? value is the probability value
vector<double>data_dedx;
TH1D *data = new TH1D("data","data",2000,0,2000);

void setdedx_distfile(const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=-99.,d_dedx=-99.;
  //  string rnd_num;
  std::string index,line;

  if (file.is_open()){
    getline(file,line); //skip header
    while(getline (file,line) ){
      
      std::istringstream in(line);
      in>>d_dedx;
      in>>d_value;
      t_dedx.push_back(d_dedx);
      t_value.push_back(d_value);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;
  dedx_dist.set_points(t_dedx,t_value);

const int size = t_dedx.size();
 double *x = &t_dedx[0];
 double *y = &t_value[0];
 
  spline = gsl_spline_alloc (gsl_interp_akima, size);
  gsl_spline_init (spline, x, y, size);

  return;
}

void setdedx_datafile(const std::string& filename) {
  ifstream file (filename.c_str());
  double d_dedx=-99.;
  //  string rnd_num;
  std::string index,line;

  if (file.is_open()){
    getline(file,line); //skip header
    while(getline (file,line) ){
      
      std::istringstream in(line);
      in>>d_dedx;
      //      cout<<"data is "<<d_dedx<<endl;
      data_dedx.push_back(d_dedx);
      data->Fill(d_dedx);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  return;
}

double gslSpline (double x, void * p){
  
  double f = 0.0;
  if(x<=t_dedx.back() && x>=t_dedx.front()){
    f = gsl_spline_eval(spline, x, acc);
    //    cout<<"x pos is "<<x<< "value of gsl spline "<< f<<endl;
  }
  else f = 0.0;
  
  return f;
}
    
double InterTable(double x, void * p){
  double f= -999999;
  if(x<=t_dedx.back() && x>=t_dedx.front())
    f = dedx_dist(x);
  else{
    f = 0;
    //    cout<<"ERROR: Input value is out of interpolation limits"<<endl;
  }
  
  return f;
}

// fcn passes back f = - 2*ln(L), the function to be minimized.

void chisq(int& npar, double* deriv, double& f, double par[], int flag){
  struct f_params alpha = {par[1],par[0]};
  double result, error;
  double expected = 1;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_function F;
  //  F.function = &InterTable;
    F.function = &gslSpline;
  F.params = &alpha;

  double ntot = data->Integral();
  int data_bins = data->GetNbinsX();
  double stepsize = (t_dedx.back()-t_dedx.front());
  //  cout<<"STEPSIZE IS "<<stepsize<<endl;

  double chisq = 0.0;
  for(int i =1; i<= data_bins; ++i){
    result = 0.0;
double adcmin = data->GetBinLowEdge(i);
    double binwidth = data->GetBinWidth(i);
    double adcmax = adcmin+binwidth;
    //    double xmin = (stepsize*par[1]*i)+par[0];
    //    double xmax =  stepsize*par[1]*(i+1)+par[0];
    double xmin = (adcmin*par[1]);//+par[0];
    double xmax =  adcmax*par[1];//+par[0];
    cout<<"xmin is "<<xmin<<" xmax is "<<xmax<<endl;
  gsl_integration_qags (&F, xmin, xmax, 0, 1e-3, 10000,w, &result, &error); 
  //    gsl_integration_qag (&F, xmin, xmax, 0, 1e-3, 1000,2, w, &result, &error); 
  //printf ("result          = % .18f\n", result);
  //  printf ("exact result    = % .18f\n", expected);
  //  printf ("estimated error = % .18f\n", error);
  //  printf ("actual error    = % .18f\n", result - expected);
  //    printf ("intervals       = %zu\n", w->size);

  //--------------
  //Calculation of log likelihood
  //--------------
  double v = result*ntot;            //expected counts in bin i
  double n = data->GetBinContent(i); //data values in bin i
    cout<<"V IS "<<v<<" N is "<<n<<endl;
    if(v!=0)chisq += pow(n-v,2)/v;
  else cout<<"Log input is 0"<<endl;
  }
  cout<<"----------------"<<endl<<endl;
  f = chisq;                    // factor of -2 so minuit gets the errors right
  cout<<f<<endl;
  gsl_integration_workspace_free (w);//clear up allocated space

}                         


void fcn(int& npar, double* deriv, double& f, double par[], int flag){
  struct f_params alpha = {par[1],par[0]};
  double result, error;
  double expected = 1;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_function F;
  //  F.function = &InterTable;
    F.function = &gslSpline;
  F.params = &alpha;

  double ntot = data->Integral();
  int data_bins = data->GetNbinsX();
  double stepsize = (t_dedx.back()-t_dedx.front());
  //  cout<<"STEPSIZE IS "<<stepsize<<endl;

  double lnL = 0.0;
  for(int i =1; i<= data_bins; ++i){
    result = 0.0;
double adcmin = data->GetBinLowEdge(i);
    double binwidth = data->GetBinWidth(i);
    double adcmax = adcmin+binwidth;
    //    double xmin = (stepsize*par[1]*i)+par[0];
    //    double xmax =  stepsize*par[1]*(i+1)+par[0];
    double xmin = (adcmin*par[1]);//+par[0];
    double xmax =  adcmax*par[1];//+par[0];
    cout<<"xmin is "<<xmin<<" xmax is "<<xmax<<endl;
  gsl_integration_qags (&F, xmin, xmax, 0, 1e-3, 10000,w, &result, &error); 
  //    gsl_integration_qag (&F, xmin, xmax, 0, 1e-3, 1000,2, w, &result, &error); 
  //printf ("result          = % .18f\n", result);
  //  printf ("exact result    = % .18f\n", expected);
  //  printf ("estimated error = % .18f\n", error);
  //  printf ("actual error    = % .18f\n", result - expected);
  //    printf ("intervals       = %zu\n", w->size);

  //--------------
  //Calculation of log likelihood
  //--------------
  double v = result*ntot;            //expected counts in bin i
  double n = data->GetBinContent(i); //data values in bin i
    cout<<"V IS "<<v<<" N is "<<n<<endl;
  if(v!=0)lnL += n*log(v);
  else cout<<"Log input is 0"<<endl;
  }
  cout<<"----------------"<<endl<<endl;
  f = -2.0 * lnL;                    // factor of -2 so minuit gets the errors right
  cout<<f<<endl;
  gsl_integration_workspace_free (w);//clear up allocated space

}                         


TH1D DrawExpected(double m, double b){
  int data_bins = data->GetNbinsX();
  TH1D expect_hist = TH1D("expected_hist","expected histogram",data_bins,0,1.*data_bins);

  double result, error;
  double expected = 1;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &InterTable;
  //  F.params = &alpha;

  double ntot = data->Integral();
  double stepsize = t_dedx.back()-t_dedx.front();
  cout<<"STEPSIZE IS "<<stepsize<<endl;

  for(int i =1; i<= data_bins; ++i){
    result = 0.0;
    double adcmin = data->GetBinLowEdge(i);
    double binwidth = data->GetBinWidth(i);
    double adcmax = adcmin+binwidth;
    //    double xmin = (stepsize*m*i)+b;
    //    double xmax =  stepsize*m*(i+1)+b;
    double xmin = adcmin*m+b;
    double xmax =  adcmax*m+b;

    //    cout<<"xmin and xmax are "<<xmin<<" "<<xmax<<endl;
  gsl_integration_qags (&F, xmin, xmax, 0, 1e-3, 1000,w, &result, &error); 
  //cout<<"AFTER INT result is "<<result<<endl;
  double v = result*ntot;            //expected counts in bin i
  //  cout<<" V IS "<<v<<endl;
  expect_hist.SetBinContent(i,v);
  }

  gsl_integration_workspace_free (w);//clear up allocated space

  //  cout<<"The equation is "<< "ADC = "<< stepsize*m<<"*[keV} + "<<b<<endl;
  return expect_hist;
}                         



int main(){
  //set dedx distribution file 
  /*
setdedx_distfile("./p_dedx.dat");
setdedx_datafile("./data_dedx_p_1.dat");
*/
  setdedx_distfile("./p_anode_dist.dat");
  setdedx_datafile("./data_anode14x10.dat");
  setdedx_datafile("./data_anode12x10.dat");
  
 // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit *minuit = new TMinuit(npar);
        minuit->SetFCN(chisq);//for chi squared
  //      minuit->SetFCN(fcn);//for max like

  double par[npar];               // the start valuesdata
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];

  par[0] = .6;            // a guess at the true value.
  stepSize[0] = 0.1;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";

  par[1] = .98;            // a guess at the true value.
  stepSize[1] = .01;       // take e.g. 0.1 of start value
  minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[1] = 0;
  parName[1] = "slope";

  for (int i=0; i<npar; i++){
    minuit->DefineParameter(i, parName[i].c_str(), 
			   par[i], stepSize[i], minVal[i], maxVal[i]);
  }
  
  // Do the minimization!
  
  minuit->Migrad();       // Minuit's best minimization algorithm
  double outpar[npar], err[npar];
  for (int i=0; i<npar; i++){
    minuit->GetParameter(i,outpar[i],err[i]);
  }
  
  cout << endl << endl << endl;
  cout << "*********************************" << endl;
  cout << "Using Maximum Likelihood, the best fit parameters are: " << endl;
  cout << "   v0: " << outpar[0] << " +/- " << err[0] << endl;
  cout << "    k: (" << outpar[1] << " +/- " << err[1] << ")E-23 " << endl;
  cout << endl;

    TH1D expected;
    //            expected = DrawExpected(.98,.4);
                    expected = DrawExpected(outpar[1],0);

  TCanvas *c1 = new TCanvas(1);
  c1->SetLogy();
  data->SetLineColor(2);
  data->SetLineWidth(3);
  data->GetXaxis()->SetRangeUser(0,200);
  expected.GetYaxis()->SetRangeUser(10e-6,10e3);
  data->Draw("lp");
  expected.SetLineWidth(3);
  expected.Draw("lp same");
  c1->SaveAs("scaled_dedxdist.png");

    //  TFile *f = new TFile("dedxscale.root","RECREATE");


  return 0;
}
