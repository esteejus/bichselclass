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

tk::spline dedx_dist;
vector<double>t_dedx,t_value;//dedx value in MeV??? value is the probability value
vector<double>data_dedx;
TH1D *bichsel_data = new TH1D("bichsel_data","bichsel data",2000,0,20);

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
      int bin = bichsel_data->GetXaxis()->FindBin(d_dedx);
      bichsel_data->SetBinContent(bin,d_value);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

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
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  return;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  TH1D bindata = TH1D("bindata","bindata",2000,0,20);
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = par[1]*data_dedx.at(i)+par[0];
    bindata.Fill(s_dedx);
  }
	
  double chisq = 0.0;
  int bins = bindata.GetXaxis()->GetNbins();
  if(bins!=2000)cout<<"BINS ARE NOT EQUAL"<<endl;
  
  for (int i=1; i<=bins; ++i){
    double v = bichsel_data->GetBinContent(i);
    double n = bindata.GetBinContent(i); //scaled dedx
    chisq += pow((n-v),2)/v;
  }
  
  f= chisq;
}



int main(){
  //set dedx distribution file 
  setdedx_distfile("./p_dedx.dat");
  setdedx_datafile("./data_dedx_p_1.data");

  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit *minuit = new TMinuit(npar);
  minuit->SetFCN(fcn);

  double par[npar];               // the start valuesdata
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];

  par[0] = .5;            // a guess at the true value.
  stepSize[0] = 0.1;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";

  par[1] = .042;            // a guess at the true value.
  stepSize[1] = 0.01;       // take e.g. 0.1 of start value
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

  TH1D *s_result = new TH1D("s_result","s_result",200,0,20);
  TH1D *expected = new TH1D("expected","expected",200,0,20);


  for(int i=0;i<data_dedx.size();++i){
    s_result->Fill(outpar[1]*data_dedx.at(i)+outpar[0]);
  }
  for(int i=0;i<t_dedx.size();++i){
    expected->Fill(t_dedx.at(i));
  }


  TH1D bindata = TH1D("bindata","bindata",2000,0,20);
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = outpar[1]*data_dedx.at(i)+outpar[0];
    bindata.Fill(s_dedx);
  }

  TCanvas *c1 = new TCanvas(1);
  c1->SetLogy();
  bindata.SetLineColor(2);
  bichsel_data->Draw();
  bindata.Draw("same");
  c1->SaveAs("chisq_dedxfit.png");

  return 0;
}
