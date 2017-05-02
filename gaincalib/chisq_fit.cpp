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
TH1D *bichsel_data = new TH1D("bichsel_data","bichsel data",2000,0,20);

void setdedx_distfile(const std::string& filename) {
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
      in>>d_value;
      sum +=d_value;
      int bin = bichsel_data->GetXaxis()->FindBin(d_dedx);
      bichsel_data->SetBinContent(bin,d_value);
    }
  }
  else std::cout << "Unable to open distribution data" << std::endl;
    bichsel_data->Scale(1./bichsel_data->Integral("width"));
cout<<"SUM ISSSSSSS "<<sum<<endl;
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
  else std::cout << "Unable to open TPC data" << std::endl;

  return;
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag){

  TH1D bindata = TH1D("bindata","bindata",2000,0,20);
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = par[1]*data_dedx.at(i)+par[0];
    //    cout<<"Data is "<<s_dedx<<endl;
bindata.Fill(s_dedx);
  }
    bindata.Scale(1./bindata.Integral("width"));
  
  double chisq = 0.0;
  int bins = bindata.GetXaxis()->GetNbins();
  if(bins!=2000)cout<<"BINS ARE NOT EQUAL"<<endl;
  
  for (int i=1; i<=bins; ++i){
    double v = bichsel_data->GetBinContent(i);
    double n = bindata.GetBinContent(i); //scaled dedx
    //     cout<<"bin "<<i<<" expected is "<<v<<" data is "<<n<<endl;
	if(v!=0)chisq += pow((v-n),2)/v;
  }
  //  cout<<"CHISQ VALUE IS "<<chisq<<endl;
  f= chisq;
}



int main(){
  //set dedx distribution file 
  setdedx_distfile("./p_dedx.dat");
  setdedx_datafile("./data_dedx_p.dat");

  // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit *minuit = new TMinuit(npar);
  minuit->SetFCN(fcn);

  double par[npar];               // the start valuesdata
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];


  par[0] = .45;            // a guess at the true value.
  stepSize[0] = 0.01;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";

  par[1] = .048;            // a guess at the true value.
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


  TH1D *bindata = new TH1D("bindata","bindata",2000,0,20);
  for(int i=0;i<data_dedx.size();++i){
    double s_dedx = outpar[1]*data_dedx.at(i)+outpar[0];
    //    cout<<"Data is "<<s_dedx<<endl;
    bindata->Fill(s_dedx);
  }
  bindata->Scale(1./bindata->Integral("width"));


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); //this will disable the title for all coming histogram

  TLatex latex;
  latex.SetTextSize(.025);
  
  TLegend *leg = new TLegend(.3,.3,.85,.6,"","NDC");
  bindata->SetLineWidth(2);
  bichsel_data->SetLineWidth(4);
  //  bindata->SetLineStyle(8);
  leg->AddEntry(bindata,"TPC data #chi^{2} best fit");
  leg->AddEntry(bichsel_data,"Bichsel straggling function");
  leg->SetLineColor(0);
  
  //setting of titles
  bichsel_data->GetXaxis()->SetTitle("Energy deposited [KeV/cm]");
  bichsel_data->GetXaxis()->CenterTitle();
  bichsel_data->GetXaxis()->SetLabelSize(.05);
  bichsel_data->GetXaxis()->SetTitleSize(.06);
    
  bichsel_data->GetYaxis()->SetTitle("Probability");
  bichsel_data->GetYaxis()->CenterTitle();
  bichsel_data->GetYaxis()->SetLabelSize(.05);
  bichsel_data->GetYaxis()->SetTitleSize(.06);
  //  bichsel_data->GetYaxis()->SetTitleOffset(1.04);

  TPaveText *pave = new TPaveText(.4,.65,.9,.85,"NDC NB");
  pave->AddText("Proton 903.63 MeV/c");
  pave->AddText("P-10 gas @ 1 atm");
  pave->AddText("1.2cm segment");
  pave->SetFillColor(0);
  
  TCanvas *c1 = new TCanvas(1);
  c1->SetLeftMargin(.15);
  c1->SetBottomMargin(.15);
  //  c1->SetLogy();
   bichsel_data->Draw();
   pave->Draw();
 bindata->SetLineColor(2);
 bindata->SetLineWidth(2);
  bindata->Draw("same");
  leg->Draw();
  c1->SaveAs("bichselbin.png");

  return 0;
}
