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

gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
tk::spline dedx_dist;
vector<double>t_dedx,t_value;//dedx value in MeV??? value is the probability value
vector<double>data_dedx;

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

double InterTable(double x, void * p){
  double value=-99999.;
  if(x<=t_dedx.back() && x>=t_dedx.front())
    value = dedx_dist(x);
  else{
    value = 0;
    cout<<"ERROR: Input value is out of interpolation limits"<<endl;
  }
  
  return value;
}

// fcn passes back f = - 2*ln(L), the function to be minimized.

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  double lnL = 0.0;
    for (int i=0; i<data_dedx.size(); ++i){
    double d_dedx = data_dedx.at(i);
    double s_dedx = d_dedx*par[1]; //scaled dedx
    double v = InterTable(s_dedx-par[0]);  //par[0]= offset, par[1]=slope
    //SHOULDNT YOU PUT LOG somehwere????
    //HWEREWEREWREWRWREWREWR
    if ( v > 0.0 ) {
      lnL += log(v);
    }
    else {
      //      lnL +=0;
            cout << "WARNING -- pdf is negative!!!" << endl;
    }
  }
  f = -2.0 * lnL;         // factor of -2 so minuit gets the errors right
  
}                         // end of fcndata_dedx_p_1.datab



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

  TCanvas *c1 = new TCanvas(1);
  c1->SetLogy();
  s_result->SetLineColor(2);
  //  expected->Draw();
    s_result->Draw("same");
  c1->SaveAs("scaled_dedxdist.png");

  TFile *f = new TFile("dedxscale.root","RECREATE");
  s_result->Scale(1./s_result->Integral("width"));
  s_result->Write();
  /*
  TH1D *test = new TH1D("test","test",200,0,200);
   for(int j=0; j<300;++j){  
     double lnL = 0.0;
     for (int i=0; i<data_dedx.size(); ++i){
       cout<<"==============="<<endl<<endl;
   cout<<"Data value is "<<data_dedx.at(i)<<endl;
       double scale =(j/2000.);
       double d_dedx = data_dedx.at(i);
       double s_dedx = d_dedx*scale; //scaled dedx
       double v = InterTable(s_dedx);  //par[0]= offset, par[1]=slope
       cout<<"Scaled value is "<<s_dedx<<endl;
     //SHOULDNT YOU PUT LOG somehwere????
    //HWEREWEREWREWRWREWREWR
      cout<<"V is "<<v<<endl;
      cout<<"Log of v is "<<log(v)<<endl;
  //cout<<"Input value "<<s_dedx<<endl;
       if ( v > 0.0 ) {
      lnL += log(v);
    }
    else {
      //         lnL +=0;
	 cout<<" vi is less than 0 Input value is "<<s_dedx<<endl;
      //       cout << "WARNING -- pdf is negative!!!" << endl;
    }
       cout<<j<<" VALUE IS "<<lnL<<endl;
     }
     test->SetBinContent(j,lnL);
   }
   test->Write();
  */
  return 0;
}
