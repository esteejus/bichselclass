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
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TPaveText.h"
#include <TAxis.h>
#include <TLine.h>

using namespace std;

gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
tk::spline dedx_dist;
tk::spline dist;
struct f_params {double a; tk::spline b; double c;};
TH1D *data = new TH1D("data","data",2000,0,2000);

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
  dedx_dist.set_points(x,y);

  return;
}

void setdedx_datafile(const std::string& filename) {
  ifstream file (filename.c_str());
  double d_dedx=-99.;

  std::string index,line;

  if (file.is_open()){
    getline(file,line); //skip header
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_dedx;
      data->Fill(d_dedx);
    }
  }
  else std::cout << "Unable to open TPC data" << std::endl;

  return;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  TH1D bindata = TH1D("bindata","bindata",100,0,100);

  int bins = bindata.GetXaxis()->GetNbins();
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = par[1]*data_dedx.at(i)+par[0];
    bindata.Fill(s_dedx);
  }
  cout<<"INTEGRAL IS "<<bindata.Integral()<<endl;
  double ntot = bindata.Integral();
  for( int i=1; i<=bins;++i){
    double value = bichsel_data->GetBinContent(i);
    value *= ntot;
    bichsel->SetBinContent(i,value);
  }
  
  double vtot = bichsel->Integral();
  double chisq = 0.0;
  if(bins!=5000)cout<<"BINS ARE NOT EQUAL"<<endl;
  
  for (int i=1; i<=bins; ++i){
    double v = bichsel->GetBinContent(i);
    double n = bindata.GetBinContent(i); //scaled dedx
    if(n!=0)  cout<<"V is "<<v<<" n is "<<n<<endl;
    double error = sqrt(pow(n,2)+pow(v,2));
    if(v!=0)chisq += log(v)*n;
  }
  //  cout<<"CHISQ VALUE IS "<<chisq<<endl;
  f= -2*(chisq-vtot);
  cout<<"CHISQ is "<<f<<endl;
}
*/



int main(){

  //  SetTable("./data_anode12.dat");
  setdedx_datafile("./data_anode12_x10.dat");

  /*
  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit *minuit = new TMinuit(npar);
  minuit->SetFCN(fcn);

  double par[npar];               // the start valuesdata
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];


  par[0] = 15;            // a guess at the true value.
  stepSize[0] = .5;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";

  par[1] = 1.2;            // a guess at the true value.
  stepSize[1] = 0.1;       // take e.g. 0.1 of start value
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


  TH1D *bindata = new TH1D("bindata","bindata",100,0,100);
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = outpar[1]*data_dedx.at(i)+outpar[0];
    //    cout<<"Data is "<<s_dedx<<endl;
    bindata->Fill(s_dedx);
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); //this will disable the title for all coming histogram

  TLatex latex;
  latex.SetTextSize(.025);
  
  TLegend *leg = new TLegend(.2,.2,.85,.4,"","NDC");
  bindata->SetLineWidth(2);
  bichsel->SetLineWidth(2);
  leg->AddEntry(bindata,"TPC data #chi^{2} best fit");
  leg->AddEntry(bichsel,"Bichsel straggling funciton");
  leg->SetLineColor(0);
  
  //setting of titles
  bichsel->GetXaxis()->SetTitle("Energy deposited [KeV/cm]");
  bichsel->GetXaxis()->CenterTitle();
  bichsel->GetXaxis()->SetLabelSize(.05);
  bichsel->GetXaxis()->SetTitleSize(.06);
    
  bichsel->GetYaxis()->SetTitle("Probability");
  bichsel->GetYaxis()->CenterTitle();
  bichsel->GetYaxis()->SetLabelSize(.05);
  bichsel->GetYaxis()->SetTitleSize(.06);
  bichsel->GetYaxis()->SetTitleOffset(1.04);

  TPaveText *pave = new TPaveText(.5,.7,.9,.85,"NDC NB");
  pave->AddText("3He 1795.7 MeV/c");
  pave->AddText("P-10 gas @ 1 atm");
  pave->AddText("1.2cm segment");
  pave->SetFillColor(0);
  
  TCanvas *c1 = new TCanvas(1);
  c1->SetLeftMargin(.15);
  c1->SetBottomMargin(.15);
  //  c1->SetLogy();
  bichsel->GetXaxis()->SetRangeUser(0,200);
  //  bichsel->GetYaxis()->SetRangeUser(0,5000);
   bichsel->Draw();
   pave->Draw();
 bindata->SetLineColor(2);
 //   double s = (bichsel->Integral()/bindata->Integral());
  bindata->Draw("same");
  leg->Draw();
  c1->SaveAs("bichselbin.png");

  return 0;
  */

}
