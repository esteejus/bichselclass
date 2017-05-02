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

tk::spline dedx_dist;

vector<double>t_dedx,t_value;//dedx value in MeV??? value is the probability value
vector<double>data_dedx;

//TH1D *bichsel_exp = new TH1D("bichsel_exp","bichsel data",1000,0,1000);
//TH1D *bichsel = new TH1D("bichsel_exp","bichsel data",100,0,100);

TH1D setdedx_distfile(const std::string& filename) {

TH1D bichsel_exp = TH1D("bichsel_exp","bichsel data",1000,0,1000);
  double sum =0;
 ifstream file (filename.c_str());
  double d_value=-99.,d_dedx=-99.;
  //  string rnd_num;
  std::string index,line;

  if (file.is_open()){
    getline(file,line); //skip header
    while(getline (file,line) ){
      
      std::istringstream in(line);
      in>>d_dedx;
      //      in>>d_value;
      //      sum +=d_value;
      //      int bin = bichsel_exp->GetXaxis()->FindBin(d_dedx);
      bichsel_exp.Fill(d_dedx);
    }
  }
  else std::cout << "Unable to open distribution data" << std::endl;
    bichsel_exp.Scale(1./bichsel_exp.Integral());

  return bichsel_exp;
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
  else std::cout << "Unable to open TPC data" << std::endl;

  return;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag,TH1D &bichsel_exp){
  TH1D expect = (TH1D) bichsel_exp->Clone();
  TH1D data = TH1D("data","bichsel data",1000,0,1000);
  int bins = data.GetXaxis()->GetNbins();
  int bich_bins = expect.GetXaxis()->GetNbins();
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = par[1]*data_dedx.at(i)+par[0];
    data.Fill(s_dedx);
  }
  double ntot = data.Integral();
  //  expect.Scale(ntot);
  double val =0;

for( int i=1; i<bins;++i){
    val= expect.GetBinContent(i);
    cout<<"VALUE BEFORE "<<val<<endl;
    cout<<"NTOT "<<ntot<<endl;
       if(ntot!=0)val = val*ntot;
    //    else value =0;
    cout<<val<<endl;
    expect.SetBinContent(i,val);
    }
  
  //  double vtot = bichsel->Integral();
  double chisq = 0.0;

  if(bich_bins!=bins)cout<<"BINS ARE NOT EQUAL"<<endl;
  for (int i=1; i<=bins; ++i){
    double v = expect.GetBinContent(i);
    double n = data.GetBinContent(i); //scaled dedx
    //    if(n!=0)  cout<<"V is "<<v<<" n is "<<n<<endl;
    //    double error = sqrt(pow(n,2)+pow(v,2));
    if(v!=0)chisq += pow(v-n,2)/v;
  }
  //  cout<<"CHISQ VALUE IS "<<chisq<<endl;
  f= chisq;
  cout<<"CHISQ is "<<f<<endl;
}




int main(){
  //set dedx distribution file 
  setdedx_distfile("./data_anode12.dat");
  setdedx_datafile("./data_anode12_x10.dat");

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


  TH1D *data = new TH1D("data","data",100,0,100);
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = outpar[1]*data_dedx.at(i)+outpar[0];
    //    cout<<"Data is "<<s_dedx<<endl;
    data->Fill(s_dedx);
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); //this will disable the title for all coming histogram

  TLatex latex;
  latex.SetTextSize(.025);
  /*
  TLegend *leg = new TLegend(.2,.2,.85,.4,"","NDC");
  data->SetLineWidth(2);
    expect->SetLineWidth(2);
  leg->AddEntry(data,"TPC data #chi^{2} best fit");
  leg->AddEntry(expect,"Expect straggling funciton");
  leg->SetLineColor(0);
  
  //setting of titles
  expect->GetXaxis()->SetTitle("Energy deposited [KeV/cm]");
  expect->GetXaxis()->CenterTitle();
  expect->GetXaxis()->SetLabelSize(.05);
  expect->GetXaxis()->SetTitleSize(.06);
    
  expect->GetYaxis()->SetTitle("Probability");
  expect->GetYaxis()->CenterTitle();
  expect->GetYaxis()->SetLabelSize(.05);
  expect->GetYaxis()->SetTitleSize(.06);
  expect->GetYaxis()->SetTitleOffset(1.04);

  TPaveText *pave = new TPaveText(.5,.7,.9,.85,"NDC NB");
  pave->AddText("3He 1795.7 MeV/c");
  pave->AddText("P-10 gas @ 1 atm");
  pave->AddText("1.2cm segment");
  pave->SetFillColor(0);
  
  TCanvas *c1 = new TCanvas(1);
  c1->SetLeftMargin(.15);
  c1->SetBottomMargin(.15);
  //  c1->SetLogy();
  expect->GetXaxis()->SetRangeUser(0,200);
  //  expect->GetYaxis()->SetRangeUser(0,5000);
  expect->Draw();
  pave->Draw();
  data->SetLineColor(2);
  //   double s = (expect->Integral()/data->Integral());
  //  data->Draw("same");
  leg->Draw();
  c1->SaveAs("expectbin.png");
  */
  return 0;
}
