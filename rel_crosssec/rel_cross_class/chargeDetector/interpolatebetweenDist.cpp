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
#include "TRandom3.h"
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

int main(){
  TFile *f = new TFile("analyzedtable.root");
  TTree *t = (TTree *)f->Get("tree");

  double bg, slope, offset;
  TH1D *refdist = nullptr;
  TGraph *scaledref = nullptr;
  TGraph *scalingrel = nullptr;
  TH1D *cdist = nullptr;
  
  t->SetBranchAddress("bg",&bg);
  t->SetBranchAddress("slope",&slope);
  t->SetBranchAddress("offset",&offset);
  t->SetBranchAddress("refdist",&refdist);
  t->SetBranchAddress("scaledref",&scaledref);
  t->SetBranchAddress("scalingrel",&scalingrel);
  t->SetBranchAddress("cdist",&cdist);

  double bg_interp = 0;
  double bg_low = 0, bg_high = 0;
  double slope_i, slope_l, slope_h;
  double offset_i, offset_l, offset_h;
  TH1D *actual_dist = nullptr;
    

  vector<double> bg_vec, slope_vec, offset_vec;

  for(int iEntry = 0;iEntry < t->GetEntries(); iEntry++)
    {
      if(iEntry%2==0)
	{
	  cout<<iEntry<<endl;
	  t->GetEntry(iEntry);
	  bg_vec.push_back(bg);
	  slope_vec.push_back(slope);
	  offset_vec.push_back(offset);
	}
    }

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *slope_table = gsl_spline_alloc(gsl_interp_linear,slope_vec.size());
  gsl_spline_init(slope_table,bg_vec.data(),slope_vec.data(),slope_vec.size());

  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  gsl_spline *offset_table = gsl_spline_alloc(gsl_interp_linear,offset_vec.size());
  gsl_spline_init(offset_table,bg_vec.data(),offset_vec.data(),offset_vec.size());


  for(int iEntry = 0;iEntry < t->GetEntries(); iEntry++)
    {
      if(!(iEntry%2==0))
	{
	  cout<<"Inside interpolation "<<endl;
	  cout<<iEntry<<endl;
	  t->GetEntry(iEntry);
	  double slope_est = gsl_spline_eval(slope_table,bg,acc);
	  double offset_est = gsl_spline_eval(offset_table,bg,acc1);
	  cout<<bg<<" "<<slope<<" "<<offset<<" "<<slope_est<<" "<<offset_est<<endl;
	  TGraph *predicted = ScaleHist(refdist,slope_est,offset_est);
	  TCanvas *c1 = new TCanvas("c1","c1",1);
	  c1->SetLogy();
	  c1->SetLogx();
	  cdist->Draw();
	  predicted->SetMarkerStyle(20);
	  predicted->SetMarkerSize(.5);
	  predicted->Draw("same PO");
	  c1->SaveAs(Form("test_interp_linear%i.png",iEntry));
	  
	}

    }

  gsl_interp_accel_free(acc);
  gsl_spline_free(slope_table);
  gsl_interp_accel_free(acc1);
  gsl_spline_free(offset_table);

  return 0;
}
