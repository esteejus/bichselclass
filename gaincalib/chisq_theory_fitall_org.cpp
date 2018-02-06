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


//struct comb_vec { vector<double> x; vector<double> y; };
struct combineData {
  vector<vector<double>> data;     
  vector<vector<double>> theory;   
  vector<vector<double>> expected; 
};

vector<struct combineData> global_container;
vector<double> par_slope, par_b, fcn_value;//global_v is the data measured
//global_f_x is expected bin content value of our model
//global vectors are the appended vectors of many theories. elment by element they are comparable


//vector<double> global_expected, global_data;

void SetTheoryVec (const std::string& filename, vector<vector<double>> &theory){

  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      theory.at(0).push_back(d_energy);
      theory.at(1).push_back(d_value);
    }
  }
  else std::cout << "Unable to open Straggling theory "<< filename << std::endl;

  return;
}

void Hist2DataVec (TH1D *hist, vector<vector<double>> &data){

  int numbins = hist->GetXaxis()->GetNbins();
  for(int i = 1; i<= numbins; i++)
    {
      double mean = hist->GetBinCenter(i);
      double value = hist->GetBinContent(i);
      double error_y = hist->GetBinError(i);
      data.at(0).push_back(mean);
      data.at(1).push_back(value);
    }


  return;
}

void SetDataVec (const std::string& filename, vector<vector<double>> &data){
  data.clear();
  data.clear();

  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;

      data.at(0).push_back(d_energy);
      data.at(1).push_back(d_value);
    }
  }
  else std::cout << "Unable to open data file" << std::endl;

  return;
}

double InterpolateTheory(double x, void *p)
{
  int idx  = *(double *) p;//idx of the global_container to get the theory from

  vector<double> theory_x = global_container.at(idx).theory.at(0);//get the x values of theory;
  vector<double> theory_y = global_container.at(idx).theory.at(1);//get the y values of theory;
  
  int size_x = theory_x.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,theory_x.data(),theory_y.data(),size_x);

  double value = 0.;
  if(x >= theory_x.front() && x <= theory_x.back())
       value = gsl_spline_eval(dist_table,x,acc);
  
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);
      
  return value;
}

void GetExpectedVec(double par[])
{

  for( int idx = 0; idx < global_container.size(); idx++)
    {
      vector<double> data_x = global_container.at(idx).data.at(0);//get the x values of theory;
      vector<double> data_y = global_container.at(idx).data.at(1);//get the y values of theory;

      vector<vector<double>>expected(2,vector<double>(0));
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
      
      double norm = accumulate(data_y.begin(), data_y.end(), 0.0);
      int ntot = data_x.size(); //number of even spaced hist bins
      double x_l_exp = data_x.front();
      double x_h_exp = data_x.back();
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
	  int alpha = idx;
	  F.function = &InterpolateTheory;
	  F.params = &alpha;
	  gsl_set_error_handler_off();
	  
	  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	  int status = gsl_integration_qag (&F, x_low, x_high, 0, 1e-3, 10000,3,w, &result, &error);    
	  if(status != 0)
	    result = 1e-10;
	  //      cout<<iBin<<endl;
	  //      cout<<"status is "<<status<<endl;
	  expected.at(0).push_back(x_center);
	  expected.at(1).push_back(result*norm);
	  gsl_integration_workspace_free (w);
	  gsl_interp_accel_free(acc);
	}
      global_container.at(idx).expected = expected;
    }//idx loop
  
  return;
}

void chisq(int& npar, double* deriv, double& f, double par[], int flag){
  //global_data
  //global_expected = GetExpected(par, data_x,data_y,theory_x,theory_y,
  GetExpectedVec(par);

  double chisq = 0.;
  
  for(int idx = 0; idx < global_container.size(); idx++)
    {
      int data_size = global_container.at(idx).data.at(1).size();
      int expected_size = global_container.at(idx).expected.at(1).size();
      if(data_size != expected_size)
	cout<<"Size of vectors not same!!!!"<<endl;

      for(int i = 0; i < data_size; i++)
	{
	  double expected_y = global_container.at(idx).expected.at(1).at(i);
	  double data_y = global_container.at(idx).data.at(1).at(i);
	  double error = sqrt(pow(data_y,2) + pow(expected_y,2));
	  //	  if(error> 0)
	    chisq += pow(expected_y - data_y, 2);

	}
    }

  //  par_b.push_back(par[0]);
  //  par_slope.push_back(par[1]);
  //  fcn_value.push_back(1./chisq);

  f = chisq;

  return;  
}

void maxll(int& npar, double* deriv, double& f, double par[], int flag){
  //global_data
  //global_expected = GetExpected(par, data_x,data_y,theory_x,theory_y,
  GetExpectedVec(par);

  double logl = 0.;
  
  for(int idx = 0; idx < global_container.size(); idx++)
    {
      int data_size = global_container.at(idx).data.at(1).size();
      int expected_size = global_container.at(idx).expected.at(1).size();
      if(data_size != expected_size)
	cout<<"Size of vectors not same!!!!"<<endl;

      for(int i = 0; i < data_size; i++)
	{
	  double expected_y = global_container.at(idx).expected.at(1).at(i);
	  double data_y = global_container.at(idx).data.at(1).at(i);

	  if(expected_y < 1.e-22)
	    {
	      	      logl += data_y*-5000.;
	      //	      cout<<"expected value is too low check fcn value"<<endl;
	    }
	  else
	    logl += data_y*log(expected_y);
	}
    }
      f = -2*logl;

  par_b.push_back(par[0]);
  par_slope.push_back(par[1]);
  fcn_value.push_back(logl);


  return;  
}


TGraph* plotTheory(int idx,double b, double slope)
{
  vector<double>theory_x,theory_y;
  theory_x = global_container.at(idx).theory.at(0);
  theory_y = global_container.at(idx).theory.at(1);

  int npoints = theory_x.size();
 TGraph *graph_temp = new TGraph(npoints);
 TGraph *graph = new TGraph(npoints);
 for(int i = 0; i<theory_x.size();i++)
   {
     double x_value = theory_x.at(i)*slope + b;
     double value = theory_y.at(i);
     graph_temp -> SetPoint(i,x_value,value);
   }

 double norm = graph_temp->Integral(1,theory_x.size());
 //calculate norm after transformation of x axis

 for(int i = 0; i<theory_x.size();i++)
   {
     double x_value = theory_x.at(i)*slope + b;
     double value = theory_y.at(i)/norm;//need to renormalize 
     graph -> SetPoint(i,x_value,value);
   }

 return graph;
}


TGraphErrors* plotData(int idx)
{
  vector<double>data_x,data_y;
  data_x = global_container.at(idx).data.at(0);
  data_y = global_container.at(idx).data.at(1);

  TGraphErrors *graph = new TGraphErrors(data_x.size());//,data_x.data(),data_y.data());
 for(int i = 0; i<data_x.size();i++)
   {
     graph->SetPoint(i,data_x.at(i),data_y.at(i));
     //     graph->SetPointError(i,0,data_error.at(i));
     //     cout<<data_y.at(i)<<endl;
   }

  return graph;
}

int main()
{


  vector<vector<double>> theory1(2,vector<double>(0));
  vector<vector<double>> data1(2,vector<double>(0));
  //    vector<vector<double>> data1(2,vector<double>(0));
    
  TFile *f = new TFile("./cocktailrootfiles/pid_2dcocktail_9_8_desat_108ndf_1641_1778mev_t.root");
  TH1D *data1_hist = (TH1D *)f->Get("c_strag");
  Hist2DataVec(data1_hist,data1);
  SetTheoryVec("cdist_p10_full_t_1612_108.data",theory1);
  
  global_container.clear();

  struct combineData entry;
  entry.data   = data1;
  entry.theory = theory1;
  //do i need expected?
  global_container.push_back(entry);




  vector<vector<double>> theory2(2,vector<double>(0));
  vector<vector<double>> data2(2,vector<double>(0));
  
  TFile *g = new TFile("./cocktailrootfiles/pid_2dcocktail_9_8_desat_108ndf_1632_1768mev_d.root");
  TH1D *data2_hist = (TH1D *)g->Get("c_strag");
  Hist2DataVec(data2_hist,data2);
  SetTheoryVec("cdist_p10_full_d_1621_108.data",theory2);

  struct combineData entry2;
  entry2.data   = data2;
  entry2.theory = theory2;
  //do i need expected?
  global_container.push_back(entry2);



  vector<vector<double>> theory3(2,vector<double>(0));
  vector<vector<double>> data3(2,vector<double>(0));
  
  TFile *h = new TFile("./cocktailrootfiles/pid_2dcocktail_9_8_desat_108ndf_897_972mev_p.root");
  TH1D *data3_hist = (TH1D *)h->Get("c_strag");
  Hist2DataVec(data3_hist,data3);
  SetTheoryVec("cdist_p10_full_p_903_108.data",theory3);

  struct combineData entry3;
  entry3.data   = data3;
  entry3.theory = theory3;
  //do i need expected?
  global_container.push_back(entry3);


   vector<vector<double>> theory4(2,vector<double>(0));
  vector<vector<double>> data4(2,vector<double>(0));
  
  TFile *i = new TFile("./cocktailrootfiles/pid_2dcocktail_9_8_desat_108ndf_878_951mev_d.root");
  TH1D *data4_hist = (TH1D *)i->Get("c_strag");
  Hist2DataVec(data4_hist,data4);
  SetTheoryVec("cdist_p10_full_d_898_108.data",theory4);

  struct combineData entry4;
  entry4.data   = data4;
  entry4.theory = theory4;
  //do i need expected?
  global_container.push_back(entry4);




 // Initialize minuit, set initial values etc. of parameters.
  const int npar = 2;              // the number of parameters
  TMinuit *minuit = new TMinuit(npar);
  minuit->SetFCN(chisq);

  double par[npar];               // the start valuesdata
  double stepSize[npar];          // step sizes 
  double minVal[npar];            // minimum bound on parameter 
  double maxVal[npar];            // maximum bound on parameter
  string parName[npar];


  par[0] = 0;            // a guess at the true value.
  stepSize[0] = 1;       // take e.g. 0.1 of start value
  minVal[0] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[0] = 0;
  parName[0] = "const";

  par[1] = 15.;            // a guess at the true value.
  stepSize[1] = .1;       // take e.g. 0.1 of start value
  minVal[1] = 0;   // if min and max values = 0, parameter is unbounded.  Only set bounds if you really think it's right!
  maxVal[1] = 0;
  parName[1] = "slope";

  for (int i=0; i<npar; i++){
    minuit->DefineParameter(i, parName[i].c_str(), 
			    par[i], stepSize[i], minVal[i], maxVal[i]);
  }
  
  // Do the minimization!


  //  minuit->Migrad();       // Minuit's best minimization algorithm
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

  /*
  TH2D *par2d = new TH2D("par2d","par2d",1000, -10, 0,1000,5, 30);
  for(int i = 0 ;i < par_slope.size(); i++)
    {
      int binx = par2d->GetXaxis()->FindBin(par_b.at(i));
      int biny = par2d->GetXaxis()->FindBin(par_slope.at(i));
      double value = fcn_value.at(i);
      //      cout<<par_b.at(i)<<endl;
      par2d->SetBinContent(binx,biny,value);
    }
  */
  TCanvas *c1 = new TCanvas("c1","c1",1);
  //    par2d->Draw("colz");
  c1->SaveAs("par2d.png");
  

  double plot_b = -8.34;
  double plot_slope = 20.2;
  
  //  double plot_b = 0;
  //  double plot_slope = 17;

  TGraphErrors *data_0 = plotData(0);
  TGraph *expected_0 = plotTheory(0,plot_b-2., plot_slope);

  TGraphErrors *data_1 = plotData(1);
  TGraph *expected_1 = plotTheory(1,plot_b,plot_slope);

  TGraphErrors *data_2 = plotData(2);
  TGraph *expected_2 = plotTheory(2,plot_b+.5,plot_slope);

  TGraphErrors *data_3 = plotData(3);
  TGraph *expected_3 = plotTheory(3,plot_b,plot_slope);


  data_1->SetMarkerSize(1);
  data_1->SetMarkerStyle(20);
  expected_1->SetMarkerStyle(20);
  expected_1->SetMarkerSize(0);
  data_1->SetMarkerColor(1);
  expected_1->SetMarkerColor(2);

  TCanvas *c2 = new TCanvas("c2","c2",1);
  //  data->GetXaxis()->SetRangeUser(8,10);
  //   c2->SetLogy();
  //  c2->SetLogx();
  data_1->GetXaxis()->SetRangeUser(50,90);
  expected_1->SetLineColor(2);
  expected_1->SetLineWidth(2);
  data_1->Draw("APO");
  expected_1->Draw("same LO");

  c2->SaveAs("fittodata_data1.png");


  data_0->SetMarkerSize(1);
  data_0->SetMarkerStyle(20);
  expected_0->SetMarkerStyle(20);
  expected_0->SetMarkerSize(0);
  data_0->SetMarkerColor(1);
  expected_0->SetMarkerColor(2);

  TCanvas *c3 = new TCanvas("c3","c3",1);
  //  data->GetXaxis()->SetRangeUser(8,10);
  //   c3->SetLogy();
  //  c3->SetLogx();
  data_0->GetXaxis()->SetRangeUser(10,100);
  expected_0->SetLineColor(2);
  expected_0->SetLineWidth(2);
  data_0->Draw("APO");
  expected_0->Draw("same LO");

  c3->SaveAs("fittodata_data0.png");


  data_2->SetMarkerSize(1);
  data_2->SetMarkerStyle(20);
  expected_2->SetMarkerStyle(20);
  expected_2->SetMarkerSize(0);
  data_2->SetMarkerColor(1);
  expected_2->SetMarkerColor(2);

  TCanvas *c4 = new TCanvas("c4","c4",1);
  //  data->GetXaxis()->SetRangeUser(8,10);
  //   c4->SetLogy();
  //  c4->SetLogx();
  data_2->GetXaxis()->SetRangeUser(10,100);
   data_2->GetXaxis()->SetRangeUser(25,45);
  data_2->GetYaxis()->SetRangeUser(0,.24);
  expected_2->SetLineColor(2);
  expected_2->SetLineWidth(2);
  data_2->Draw("APO");
  expected_2->Draw("same LO");

  c4->SaveAs("fittodata_data2_p.png");

  data_3->SetMarkerSize(1);
  data_3->SetMarkerStyle(20);
  expected_3->SetMarkerStyle(20);
  expected_3->SetMarkerSize(0);
  data_3->SetMarkerColor(1);
  expected_3->SetMarkerColor(2);

  TCanvas *c5 = new TCanvas("c5","c5",1);
  //  data->GetXaxis()->SetRangeUser(8,10);
  //   c5->SetLogy();
  //  c5->SetLogx();
  data_3->GetXaxis()->SetRangeUser(70,200);
  data4_hist->GetXaxis()->SetRangeUser(85,125);
  data4_hist->GetYaxis()->SetRangeUser(0,.12);
  expected_3->SetLineColor(2);
  expected_3->SetLineWidth(2);
  data_3->Draw("APO");
  expected_3->Draw("same LO");

  c5->SaveAs("fittodata_data3_d.png");


  return 0;

}
