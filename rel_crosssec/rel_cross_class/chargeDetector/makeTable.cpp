#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <gsl/gsl_spline.h>

using namespace std;


TGraph * ScaleHist(TH1 *h, double scale, double offset)
{
  vector<double> x_vec,y_vec;
  for(int i = 0; i < h->GetNbinsX(); i++)
    {
      double x = h->GetBinCenter(i);
      double y = h->GetBinContent(i);
      x_vec.push_back(x*scale + offset);
      y_vec.push_back(y);
    }

  TGraph *temp = new TGraph(x_vec.size(),x_vec.data(),y_vec.data());
  double norm = temp->Integral(0);
    for(int i = 0; i < y_vec.size(); i++)
    y_vec.at(i) /= norm;
  TGraph *graph = new TGraph(x_vec.size(),x_vec.data(),y_vec.data());

  return graph;
}


TGraph * GetCScaling(TH1 *h1, TH1 *h2){

  vector<double> h1_x, h1_y, h2_x, h2_y;
  double h1_scale = h1->Integral("width");
  double h2_scale = h2->Integral("width");
  cout<<h1_scale<<endl;
  h1->Scale(1./h1_scale);
  h2->Scale(1./h2_scale);
  
  int binh1 = h1->GetXaxis()->GetNbins();
  int binh2 = h2->GetXaxis()->GetNbins();  

  int h1_low = h1->FindFirstBinAbove(1e-8);
  int h1_high = h1->FindLastBinAbove(1e-8);
  int h2_low = h2->FindFirstBinAbove(1e-8);
  int h2_high = h2->FindLastBinAbove(1e-8);

  double prev_cdf = -1;
  for(int i = h1_low-1 ;i < h1_high+1; i++)
    {
      double h1_x_value = h1 -> GetBinCenter(i);
      double h1_cdf = h1 -> Integral(1,i,"width");
      if(h1_cdf > prev_cdf)//this ensures monotomic 
	{
	  h1_x.push_back(h1_x_value);
	  h1_y.push_back(h1_cdf);
	  prev_cdf = h1_cdf;
	}
    }
  prev_cdf = -1;
  for(int i = h2_low-1 ;i < h2_high+1; i++)
    {
      double h2_x_value = h2 -> GetBinCenter(i);
      double h2_cdf = h2 -> Integral(1,i,"width");
      if(h2_cdf > prev_cdf)//this ensures monotomic 
	{
	  h2_x.push_back(h2_x_value);
	  h2_y.push_back(h2_cdf);
	  prev_cdf = h2_cdf;
	}
    }

int size_h2 = h2_y.size();
gsl_interp_accel *acc_h1 = gsl_interp_accel_alloc();
gsl_spline *h2_table = gsl_spline_alloc(gsl_interp_akima,size_h2);
gsl_spline_init(h2_table,h2_y.data(),h2_x.data(),size_h2);

int size_h1 = h1_y.size();
gsl_interp_accel *acc_h2 = gsl_interp_accel_alloc();
gsl_spline *h1_table = gsl_spline_alloc(gsl_interp_akima,size_h1);
gsl_spline_init(h1_table,h1_y.data(),h1_x.data(),size_h1);
 

 int npoints = 1e3;
 TGraph *graph = new TGraph(npoints);
 double step_size = 1./npoints;
 for(int i = 0 ;i < npoints; i++)
   {
     double x = step_size*i;
     double energy_h1 = gsl_spline_eval(h1_table,x,acc_h1);
     double energy_h2 = gsl_spline_eval(h2_table,x,acc_h2);
     graph -> SetPoint(i,energy_h1,energy_h2);
   }

 //  TGraph *graph = new TGraph(h1_x.size(),h1_x.data(),h1_y.data());

  return graph;
}

int main(){

  TFile *input = new TFile("inputtable.root");
  TTree *tree = (TTree *)input->Get("table");
  double bg = 0;
  int ns = 0;
  TH1D *hist = nullptr;

  auto b_bg  = tree->GetBranch("bg");
  auto b_cdist  = tree->GetBranch("cdist");
  b_bg->SetAddress(&bg);
  b_cdist->SetAddress(&hist);

  const int num_entries = tree->GetEntries();

  double bg_ary[num_entries];
  TH1D *cdist[num_entries];

  for(int iEntry = 0; iEntry < num_entries; iEntry++)
    {
      tree->GetEntry(iEntry);
      cdist[iEntry]  = (TH1D*)hist->Clone();
      bg_ary[iEntry] = bg;
    }

  TH1D *refdist = cdist[num_entries/2];//SET REFERENCE FUNCTION HERE!
  TGraph *scaling[num_entries];
  TGraph *c_scaling = nullptr;
  TGraph *scaledref = nullptr;
  TH1D *c_cdist     = nullptr;
  
  TF1 *line = new TF1("line","[0]*x + [1]",0,1e5);
  line->SetParameters(-100,10);

  //      scaling[0] = GetCScaling(refdist,refdist);
  //      scaling[0]->Fit("line");

  double slope = 0, offset = 0;

  TFile *out = new TFile("analyzedtable.root","RECREATE");
  out->cd();
  TTree *t = new TTree("tree","tree");
  t->Branch("bg",&bg);
  t->Branch("slope",&slope);
  t->Branch("offset",&offset);
  t->Branch("refdist",&refdist);
  t->Branch("cdist",&c_cdist);
  t->Branch("scaledref",&scaledref);
  t->Branch("scalingrel",&c_scaling);


  for(int iEntry = 0; iEntry < num_entries; iEntry++)
    {
      scaling[iEntry] = GetCScaling(refdist,cdist[iEntry]);
      scaling[iEntry]->Fit("line");
      slope = line->GetParameter(0);
      offset = line->GetParameter(1);
      TGraph *graph = ScaleHist(refdist,slope,offset);

      bg = bg_ary[iEntry];
      c_cdist = (TH1D *)cdist[iEntry]->Clone();
      scaledref = (TGraph *)graph->Clone();
      c_scaling = (TGraph *)scaling[iEntry]->Clone();
      t->Fill();
  
      /*
      TCanvas *c1 = new TCanvas("c1","c1",1);
      c1->SetLogy();
      c1->SetLogx();
      //      cdist[iEntry]->GetXaxis()->SetRangeUser(0,1e4);
      cdist[iEntry]->Draw();
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(.5);
      graph->Draw("same PO");
      //      scaling[iEntry]->SetMarkerStyle(22);
      //      scaling[iEntry]->SetMarkerSize(2);
      //      scaling[iEntry]->Draw("APO");
      c1->SaveAs(Form("%i.png",iEntry));
      */
    }  
  t->Write();
  out->Close();

  return 0;

}
