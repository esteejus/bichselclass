#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <TMinuit.h>

using namespace std;

vector<double> strag_vec,strag_value;
vector<double> data_vec,data_value,data_error;
vector<double> expected_vec,expected_value;

void SetTheoryVec (const std::string& filename) {
  strag_vec.clear();
  strag_value.clear();
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      strag_vec.push_back(d_energy);
      strag_value.push_back(d_value);
    }
  }
  else std::cout << "Unable to open Straggling theory" << std::endl;

  return;
}

void Hist2DataVec (TH1D *hist) {
  data_vec.clear();
  data_value.clear();
  int numbins = hist->GetXaxis()->GetNbins();
  for(int i = 1; i<= numbins; i++)
    {
      double mean = hist->GetBinCenter(i);
      double value = hist->GetBinContent(i);
      double error_y = hist->GetBinError(i);
      data_vec.push_back(mean);
      data_value.push_back(value);
      data_error.push_back(error_y);
    }

  return;
}

void SetDataVec (const std::string& filename) {
  data_vec.clear();
  data_value.clear();

  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;

      data_vec.push_back(d_energy);
      data_value.push_back(d_value);
    }
  }
  else std::cout << "Unable to open data file" << std::endl;

  return;
}

double InterpolateTheory(double x, void *)
{
  int size_x = strag_vec.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,strag_vec.data(),strag_value.data(),size_x);

  double value = 0.;
  if(x >= strag_vec.front() && x <= strag_vec.back())
       value = gsl_spline_eval(dist_table,x,acc);
  
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);
      
  return value;
}

void GetExpectedVec(double par[])
{

  expected_vec.clear();
  expected_value.clear();
  //In data with ntot events
  // x_l is lower bin edge
  // x_h is uppber bin edge
  // with bin width w = (x_h - x_l)/ntot
  //Relation to theory is x_exp = x_th *slope + b
  //In the theory then
  // x_l = x_l' * slope + b
  // h = (x_h - x_l)/ntot
  // h * ntot = (x_h' *slope + b) - (x_l'*slope +b)
  // h * (ntot/slope) = x_h' - x_l'
  // x_h' =  h *( ntot/slope) + x_l'

  double norm = accumulate(data_value.begin(), data_value.end(), 0.0);
  int ntot = data_vec.size(); //number of even spaced hist bins
  double x_l_exp = data_vec.front();
  double x_h_exp = data_vec.back();
  double h_exp = (x_h_exp - x_l_exp)/(ntot-1);
  //x_l_exp and x_h_exp are centers of histogram bins
  //thus distance between two will be (ntot-1)*h_exp
  x_l_exp -= (h_exp/2);//subtract half a bin to get low edge
  x_h_exp += (h_exp/2);//add half bin to get upper edge
  
  double x_l_th = (x_l_exp - par[0])/par[1];
  double x_h_th = h_exp * (ntot/par[1]) + x_l_th;
  double h_th = (x_h_th - x_l_th)/ntot;
  //here x_l_th and x_h_th are correct low and high edge
  //thus the distance is ntot bins
  
  for(int iBin = 0 ;iBin < ntot; iBin++)
    {

      gsl_interp_accel *acc = gsl_interp_accel_alloc();
      double x_low = h_th * iBin + x_l_th;
      double x_high = h_th * (iBin + 1) + x_l_th;
      double x_center = x_low + (h_th/2.);
      double result = 0., error = 0.;
      gsl_function F;
      F.function = &InterpolateTheory;
      gsl_set_error_handler_off();

      gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
      int status = gsl_integration_qag (&F, x_low, x_high, 0, 1e-3, 10000,3,w, &result, &error);    
      if(status != 0)
	result = 1e-10;
      //      cout<<iBin<<endl;
      //      cout<<"status is "<<status<<endl;
      expected_vec.push_back(x_center);
      expected_value.push_back(result*norm);

      gsl_integration_workspace_free (w);
      gsl_interp_accel_free(acc);
    }

  return;
}

/*
void fcn(int& npar, double* deriv, double& f, double par[], int flag){
  GetExpectedVec(par);
  double chisq = 0.;
  if(expected_vec.size() != data_vec.size())
    cout<<"Size of vectors not same!!!!"<<expected_vec.size()<<endl;
  
    for(int i = 0; i < expected_vec.size(); i++)
      {
	if(expected_value.at(i) < 1e-10)
	  {
	    //chisq += 1e10;
	    cout<<"expected value is too low check fcn value"<<endl;
	  }
	else
	  chisq += pow(data_value.at(i)-expected_value.at(i),2)/expected_value.at(i);
      }
    f = chisq;

  return;  
}
*/



void fcn(int& npar, double* deriv, double& f, double par[], int flag){
  GetExpectedVec(par);
  double logl = 0.;
  if(expected_vec.size() != data_vec.size())
    cout<<"Size of vectors not same!!!!"<<expected_vec.size()<<endl;
  
    for(int i = 0; i < expected_vec.size(); i++)
      {
	if(expected_value.at(i) < 1.e-20)
	  {
	    logl += 0;
	    cout<<"expected value is too low check fcn value"<<endl;
	  }
	else
	  logl += data_value.at(i)*log(expected_value.at(i));
      }
    f = -2*logl;
    //    cout<<"f "<<f<<endl;
  return;  
}


TGraph* plotTheory(double b, double slope)
{
  int npoints = strag_vec.size();
 TGraph *graph_temp = new TGraph(npoints);
 TGraph *graph = new TGraph(npoints);
 for(int i = 0; i<strag_vec.size();i++)
   {
     double x_value = strag_vec.at(i)*slope + b;
     double value = strag_value.at(i);
     graph_temp -> SetPoint(i,x_value,value);
   }

 double norm = graph_temp->Integral(1,strag_vec.size());
 //calculate norm after transformation of x axis

 for(int i = 0; i<strag_vec.size();i++)
   {
     double x_value = strag_vec.at(i)*slope + b;
     double value = strag_value.at(i)/norm;//need to renormalize 
     graph -> SetPoint(i,x_value,value);
   }

 return graph;
}


TGraphErrors* plotData()
{
  TGraphErrors *graph = new TGraphErrors(data_vec.size());//,data_vec.data(),data_value.data());
 for(int i = 0; i<data_vec.size();i++)
   {
     graph->SetPoint(i,data_vec.at(i),data_value.at(i));
     graph->SetPointError(i,0,data_error.at(i));
     //     cout<<data_value.at(i)<<endl;
   }

  return graph;
}

int main()
{

  //  SetDataVec("test_data.dat");
  //  SetTheoryVec("test_theory.dat");

  TFile *f = new TFile("./cocktailrootfiles/pid_2dcocktail_9_8_desat_108ndf_1641_1778mev_t.root");
  //  TH1D *data_in = (TH1D *)f->Get("full_strag");
  TH1D *data_in = (TH1D *)f->Get("full_strag");
  Hist2DataVec(data_in);
  SetTheoryVec("p10_full_t_1612_108.data");

 // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit *minuit = new TMinuit(npar);
  minuit->SetFCN(fcn);

  double par[npar];               // the start valuesdata
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];


  par[0] = -5;            // a guess at the true value.
  stepSize[0] = .1;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";

  par[1] = 16.;            // a guess at the true value.
  stepSize[1] = .1;       // take e.g. 0.1 of start value
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
  cout << "    k: (" << outpar[1] << " +/- " << err[1] << endl;
  cout << endl;



  TGraphErrors *data = plotData();
  TGraph *expected = plotTheory(-13.81,17.59);
  //  TGraph *expected = plotTheory(-5.54,17.21);
  //  data->SetMarkerStyle(20);
  data->SetMarkerSize(.5);
  expected->SetMarkerStyle(20);
  expected->SetMarkerSize(0);


  data->SetMarkerColor(1);
  expected->SetMarkerColor(2);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  //  data->GetXaxis()->SetRangeUser(8,10);
  //   c1->SetLogy();
  c1->SetLogx();
  data->GetXaxis()->SetRangeUser(10,500);
  expected->SetLineColor(2);
  expected->SetLineWidth(2);
  data->Draw("APO");
  expected->Draw("same LO");

  c1->SaveAs("fittodata.png");

  return 0;

}
