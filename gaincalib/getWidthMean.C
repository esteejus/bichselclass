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

void getWidthMean(){

  //  SetDataVec("test_data.dat");
  //  SetTheoryVec("test_theory.dat");

  TFile *f = new TFile("./cocktailrootfiles/pid_2dcocktail_9_8_desat_108ndf_859_930mev_t.root");
  
  //  TH1D *data_in = (TH1D *)f->Get("full_strag");
  TH1D *data_in = (TH1D *)f->Get("c_strag");
  Hist2DataVec(data_in);
  SetTheoryVec("cdist_p10_full_t_1700.data");
  TGraph *theory = plotTheory(0,1);
  
  cout<<"Mean "<<data_in->GetMean()<< "width "<<data_in->GetStdDev()<<endl;
  theory->Fit("gaus");
  theory->Draw();

  //data_in->Fit("gaus");


  
  return 0;

}
