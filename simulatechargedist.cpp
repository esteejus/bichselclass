#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "bichsel.h"
#include "TRandom3.h"

using namespace std;

void GenerateTrack(double z,double mass, double mom, double t_length, vector<vector<double>> &padplane);

vector<vector<double>> InitPadPlane()
{
 //initialize the padplane
  int num_row = 108;
  int num_layers = 112;
  vector<vector<double>> padplane(num_layers + 1);//add 1 to so we can use 112 as a valid entry
  for(int i = 0; i <= num_layers; i++)
    padplane.at(i).resize(108+1,0.);//add 1 to so we can use 108 as valid element 

  return padplane;
}

std::tuple<int,int,double,double> GetCurrentPad(double xpos, double ypos)

{
  int row = floor(xpos/8.);
  int layer = floor(ypos/12.);

  double x_inpad = (xpos - row * 8) - 4;//in coordiantes of tpc current pad
  double y_inpad = (ypos - layer * 12) - 6;//in coordinates of tpc current pad


  //this section propogates which wire in the current pad will the charge terminate on
  //there are 3 wires over each pad organized as such
  //Each A is the anode wire
  //the second anode wire is at y=0 in pad coordinates
  //<--2mm--> A <---4mm---> A <---4mm---> A <--2mm-->
  //
  //
  //
  //-------------------------------------------------
  //This is the PAD ^^^^
  //<-------- -y ----------|---------- +y ---------->
 
  if ( y_inpad < -2 && y_inpad >=-6)
    y_inpad = -4;
  else if(y_inpad < 2 && y_inpad >=-2)
    y_inpad = 0;
  else if(y_inpad <= 6 && y_inpad >=2)
    y_inpad = 4;
  else
    cout<<"Y inpad is not within bounds for some reason"<<endl;
  
  return std::make_tuple(row,layer,x_inpad,y_inpad);
}

double GetFract(double x1, double x2, double y1, double y2, double x_0, double y_0)
{
  //x1 x2 integral in row direciton from x1 to x2 in current pad reference frame
  //y1 y2 integral in layer direciton from y1 to y2 in current pad reference frame
  double sigma_x = .1;
  double sigma_y = .1;
  //  double sigma_x = 5;
  //  double sigma_y = 5;

  double u1 = (x1 - x_0)/sigma_x;
  double u2 = (x2 - x_0)/sigma_x;
  double v1 = (y1 - y_0)/sigma_y;
  double v2 = (y2 - y_0)/sigma_y;

  return (1./4)*(erf(v2)-erf(v1))*(erf(u2)-erf(u1));
}

vector<std::tuple<int,int,double>> GetPadResponse(std::tuple<int,int,double,double> curPad , double charge)
{
  vector<std::tuple<int,int,double>> unitpad;
  double x_step = 8;
  double y_step = 12;

  int row = std::get<0>(curPad);
  int layer = std::get<1>(curPad);
  double x_0 = std::get<2>(curPad);
  double y_0 = std::get<3>(curPad);

  //making a 3x3 grid of unit pad plane
  //pad row = 0 layer = 0 is current pad
  for(int i = -2; i <= 2; i++)//row
    {
      for(int j = -1; j <= 1; j++)//layer
      	{
	  double y1 = -6 + j*y_step;
	  double y2 =  6 + j*y_step;
	  double x1 = -4 + i*x_step;
	  double x2 =  4 + i*x_step;

	  if(row + i > 0 && layer + j > 0)
	    {
	      unitpad.emplace_back(row + i,layer + j,charge*GetFract(x1,x2,y1,y2,x_0,y_0));
	      cout<<"x1 x2 y1 y2 "<<x1<<" "<<x2<<endl;
	      cout<<row + i<<" "<<layer + j<<" "<< charge*GetFract(x1,x2,y1,y2,x_0,y_0)<<endl;
	    }
	}
    }

  return unitpad;
}

void FillPadPlane(vector<std::tuple<int,int,double>> hits, vector<vector<double>> &padplane)
{
  for(int iHit = 0; iHit < hits.size(); iHit++)
    {
      int row = std::get<0>(hits.at(iHit));
      int layer = std::get<1>(hits.at(iHit));
      double charge = std::get<2>(hits.at(iHit));
      //      cout<<std::get<2>(hits.at(iHit))<<endl;
      if(charge>0)
	{
	  cout<<"Row layer "<<row<<" "<<layer<<endl;
	  cout<<"pad plane charge is "<<      padplane.at(layer).at(row)<<endl;
	  cout<<"charge is "<<charge<<endl;
	}
    padplane.at(layer).at(row) += charge;
    }
  
  return;
}

void PlotClusterDist(vector<vector<double>>padplane,TH1D *dedx_dist)
{
  //  TH1D *dedx_dist = new TH1D("clusterdist","clusterdist",1000,0,1e5);

  int num_layers = 112;
  int num_rows = 108;
  for(int iLayer = 1; iLayer<= num_layers; iLayer++)
    {
      vector<double> curLayer = padplane.at(iLayer);
      curLayer.erase(remove_if(curLayer.begin(),curLayer.end(),[](double x){ return x < 1e-2;}),curLayer.end());
      double sum = std::accumulate(curLayer.begin(),curLayer.end(),0.0);
      if(sum > 0)
	{
	  //	  cout<<sum<<endl;
	  dedx_dist -> Fill(sum);
	}
    }

  return;
}

double GetC(vector<vector<double>>padplane, double threshold, double trunc_f)
{
  int num_layers = 112;
  int num_rows = 108;
  vector<double> chg_layers;
  for(int iLayer = 1; iLayer<= num_layers; iLayer++)
    {
      vector<double> curLayer = padplane.at(iLayer);
      curLayer.erase(remove_if(curLayer.begin(),curLayer.end(),[threshold](double x){ return x<threshold;}),curLayer.end());
      double sum = std::accumulate(curLayer.begin(),curLayer.end(),0.0);
      if(sum > 0)
	chg_layers.push_back(sum/1.2);
    }

  std::sort(chg_layers.begin(),chg_layers.end(),[](double i, double j){return i < j;});
  int element = trunc_f*chg_layers.size();
  element -= 1;//erase saves the elment you put in. to erase element a this pos sub 1
  chg_layers.erase(chg_layers.begin()+element,chg_layers.end());
  double sum = std::accumulate(chg_layers.begin(),chg_layers.end(),0.0);
  sum = sum/chg_layers.size();
  
  return sum;
}

void GetCDist(double z, double mass, double mom, double t_length, double threshold, double trunc_f, int n_events, TH1D * c_dist)
{
  //void GenerateTrack(double z,double mass, double mom, double t_length, vector<vector<double>> &padplane)
  //double GetC(vector<vector<double>>padplane, double threshold, double trunc_f)
  //  TH1D *c_dist = new TH1D("c_dist","c_dist",100,0,1e5);
  for(int iEvent = 0 ;iEvent < n_events; iEvent++)
    {
      if(iEvent%10==0)
	cout<<"Event number: "<<iEvent<<endl;
      vector<vector<double>> padplane = InitPadPlane();
      GenerateTrack(z,mass,mom,t_length,padplane);
      double c_avg = GetC(padplane,threshold,trunc_f);
      c_dist -> Fill(c_avg);
    }
  
  return;
}

TH2D *PlotPadPlane(vector<vector<double>> padplane)
{
  TH2D *padplane_hist = new TH2D("","",112,0,112,108,0,108);
  int num_layers = 112;
  int num_rows = 108;
  for(int iLayer = 1; iLayer<= num_layers; iLayer++)
    {
      vector<double> curLayer = padplane.at(iLayer);
	for(int iRow = 1;iRow <= num_rows;iRow++)
	{
	  double charge = curLayer.at(iRow);
	  padplane_hist->SetBinContent(iLayer,iRow,charge);
	  //	  cout<<charge<<endl;
	}
    }

  return padplane_hist;
}

void Plot2DPID(double p_low, double p_high, double z,double mass, double t_length,double threshold,double trunc_f, int num_events,TH2D *pid)
{
  TRandom3 *ran = new TRandom3(0);
  for(int iEvent = 0; iEvent < num_events; iEvent++)
    {
      double mom = ran->Uniform(p_low,p_high);
      vector<vector<double>> padplane = InitPadPlane();
      GenerateTrack(z,mass,mom,t_length,padplane);
      double c_avg = GetC(padplane,threshold,trunc_f);
      //      cout<<mom<<" "<<c_avg<<endl;
      pid->Fill(mom,c_avg);
    }

  return;
}

void GenerateTrack(double z,double mass, double mom, double t_length, vector<vector<double>> &padplane)
{
  double b_gamma = mom/mass;
  Bichsel b(b_gamma);
  b.SetInvXSec("P10M0invw_31623.inv");
  b.SetM0Table("m0.dat");

  double x_0 = 432;//initial x pos approx center of pad plane
  double y_0 = 0; //initial y pos really z axis on real TPC 
  double dist = 0; //distance traveled 
  while(dist < t_length)
    {
      double dx = b.GetMCstep()*10;//GetMCStep is in [cm]
      double dE = b.GetEloss();
      dist +=  dx;
      auto hit = GetCurrentPad(432,dist);//find hit from x,y info
      cout<<"Charge of event is "<<dE*pow(z,2)<<" "<<dist<<endl;
      auto hit_response = GetPadResponse(hit,dE*pow(z,2));
      FillPadPlane(hit_response,padplane);
    }

  return;
}

Double_t median1(TH1D *h1)
{
  //compute the median for 1-d histogram h1
  Int_t nbins = h1->GetXaxis()->GetNbins();
  Double_t *x = new Double_t[nbins];
  Double_t *y = new Double_t[nbins];
  for (Int_t i=0;i<nbins;i++) {
    x[i] = h1->GetXaxis()->GetBinCenter(i+1);
    y[i] = h1->GetBinContent(i+1);
  }
  Double_t median = TMath::Median(nbins,x,y);
  delete [] x;
  delete [] y;
  return median;
}

TGraph * median2(TH2D *h2, double z, double mass)
{
  //compute and print the median for each slice along X of h2
  Int_t nbins = h2->GetXaxis()->GetNbins();
  TGraph *output = new TGraph(nbins);
  TString name = h2 -> GetName();
  name += "_mean";
  double lowbin = h2 -> GetXaxis() -> GetBinLowEdge(1);
  double highbin = h2 -> GetXaxis() -> GetBinUpEdge(nbins);

  for (Int_t i=1;i<=nbins;i++) {
    TH1D *h1 = h2->ProjectionY("i",i,i);
    double mom = h2 -> GetXaxis() -> GetBinCenter(i);
    Double_t mean = h1->GetMean();
    double error = h1->GetEntries();
    error = mean/sqrt(error);
    if(mean > 1e-1){
      output -> SetPoint(i,mom/mass,mean/pow(z,2));
      //      output -> SetBinError(i,error);
    }
    //    printf("Median of Slice %d, Median=%g, Mean = %g\n",i,median,mean);
    delete h1;
  }

  return output;
}



int main()
{
  vector<vector<double>> padplane = InitPadPlane();
 

  double z       = 1;
  double mass    = 938;
  double mom     = 2814;
  double length  = 12;//mm
  double thresh  = 1;//adc not right now tho
  double trunc_f = .7;
  int    events  = 1e3;

  GenerateTrack(z,mass,mom,length,padplane);
  TH2D *pp = PlotPadPlane(padplane);
  TH1D *dedx_dist = new TH1D("clusterdist","clusterdist",100,0,1e4);
  PlotClusterDist(padplane,dedx_dist); 

  TCanvas *c3 = new TCanvas("c3","c3",1);
  pp->GetXaxis()->SetRangeUser(0,3);
  pp->GetYaxis()->SetRangeUser(50,57);
  pp->Draw("colz");
  c3->SaveAs("padplane.png");
  //  cout<< GetC(padplane,1e-3,.7)<<endl;

  //TH1D * GetCDist(double z, double mass, double mom, double t_length, double threshold, double trunc_f, int n_events)
  //  TH1D *c_dist = new TH1D("c_dist","c_dist",200,0,5e3);
  //  GetCDist(z,mass,mom,length,thresh,trunc_f,events,c_dist);
  /*
  //  auto padplane = GetPadResponse(d,1);
  cout<<  std::get<0>(padplane.at(0)) <<endl;
  cout<<  std::get<1>(padplane.at(0)) <<endl;    
  cout<<  std::get<2>(padplane.at(0)) <<endl;

  TH2D *unitpad = new TH2D("","",5,0,5,5,0,5);
  for(int i = 0; i < padplane.size(); i++)
    {
      int row = std::get<0>(padplane.at(i));
      int layer = std::get<1>(padplane.at(i));
      double charge = std::get<2>(padplane.at(i));
      unitpad -> SetBinContent(row+3,layer+3,charge);
      cout<<row<<" "<<layer<<" "<<charge<<endl;
    }
  */

  /*
  double mom_low   = 800;
  double mom_high  = 3000;
  int num_events = 2e3;

  double t_length = 1000;
  double fract = .7;
  double p_z    = 1;
  double d_z    = 1;
  double t_z    = 1;
  double he3_z    = 2;
  double he4_z    = 2;

 double amu=931.5;     //MeV/c^2 from Nuclear Wallet cards
  //masses from nuclear wallet cards
 double p_mass = amu + 7.2889;
 double d_mass = amu*2+13.136;
 double t_mass = amu*3+14.95;
 double he3_mass = amu*3+14.931;
 double he4_mass = amu*4+2.425;
  
 
  TH2D *pid_p = new TH2D("pid","pid",100,0,5000,5000,0,20000);
  TH2D *pid_d = new TH2D("pid_d","pid_d",100,0,5000,5000,0,20000);
  TH2D *pid_t = new TH2D("pid_t","pid_t",100,0,5000,5000,0,20000);
  TH2D *pid_he3 = new TH2D("pid_he3","pid_he3",100,0,5000,5000,0,20000);
  TH2D *pid_he4 = new TH2D("pid_he4","pid_he4",100,0,5000,5000,0,20000);

  //Plot2DPID(double p_low, double p_high, double z,double mass, double t_length,double threshold,double trunc_f, int num_events,TH2D *pid)
  Plot2DPID(mom_low*p_z,mom_high*p_z,p_z,p_mass,t_length,thresh,fract,num_events,pid_p);
  cout<<"Finished proton "<<endl;
  Plot2DPID(mom_low*d_z,mom_high*d_z,d_z,d_mass,t_length,thresh,fract,num_events,pid_d);
  cout<<"Finished d "<<endl;
  Plot2DPID(mom_low*t_z,mom_high*t_z,t_z,t_mass,t_length,thresh,fract,num_events,pid_t);
  cout<<"Finished t "<<endl;
  Plot2DPID(mom_low*he3_z,mom_high*he3_z,he3_z,he3_mass,t_length,thresh,fract,num_events,pid_he3);
  cout<<"Finished he3 "<<endl;
  Plot2DPID(mom_low*he4_z,mom_high*he4_z,he4_z,he4_mass,t_length,thresh,fract,num_events,pid_he4);
  cout<<"Finished he4 "<<endl;


 //start of line
  TGraph * p_mean = median2(pid_p,p_z,p_mass);
  TGraph * d_mean = median2(pid_d,d_z,d_mass);
  TGraph * t_mean = median2(pid_t,t_z,t_mass);
  TGraph * he3_mean = median2(pid_he3,he3_z,he3_mass);
  TGraph * he4_mean = median2(pid_he4,he4_z,he4_mass);

  int p_col = kMagenta;
  int d_col = kBlue;
  int t_col = kCyan;
  int he3_col = kGreen;
  int he4_col = kRed;
  int li6_col = kOrange;
  int li7_col = kOrange+3;

  int p_mark = 20;
  int d_mark = 21;
  int t_mark = 22;
  int he3_mark = 23;
  int he4_mark = 20;
  int li6_mark = 33;
  int li7_mark = 34;

  double  p_size   = 1.75;
  double  d_size   = 1.75;
  double  t_size   = 1.75;
  double  he3_size = 1.75;
  double  he4_size = 1.75;
  double  li6_size = 2;
  double  li7_size = 2;

  p_mean -> SetMarkerColor(p_col);
  d_mean -> SetMarkerColor(d_col);
  t_mean -> SetMarkerColor(t_col);
  he4_mean -> SetMarkerColor(he4_col);
  he3_mean -> SetMarkerColor(he3_col);


  p_mean -> SetLineColor(p_col);
  d_mean -> SetLineColor(d_col);
  t_mean -> SetLineColor(t_col);
  he3_mean -> SetLineColor(he3_col);
  he4_mean -> SetLineColor(he4_col);


  p_mean -> SetMarkerStyle(p_mark);
  d_mean -> SetMarkerStyle(d_mark);
  t_mean -> SetMarkerStyle(t_mark);
  he3_mean -> SetMarkerStyle(he3_mark);
  he4_mean -> SetMarkerStyle(he4_mark);

  p_mean -> SetMarkerSize(p_size);
  d_mean -> SetMarkerSize(d_size);
  t_mean -> SetMarkerSize(t_size);
  he3_mean -> SetMarkerSize(he3_size);
  he4_mean -> SetMarkerSize(he4_size);

  TCanvas *c1 = new TCanvas("","",1);
  p_mean -> GetYaxis()->SetRangeUser(0,4000);
  p_mean -> Draw("APO");
  d_mean -> Draw("same P0");
  t_mean -> Draw("same P0");
  he3_mean -> Draw("same P0");
  he4_mean -> Draw("same P0");

  //  pp->GetYaxis()->SetRangeUser(46,60);
  //  pp ->Draw("colz");
  //  dedx_dist ->Draw();
  //  pid_p->Draw("colz");
  //  pid_d ->Draw("colz SAME");
  //  pid_t ->Draw("colz SAME");
  //  pid_he3 ->Draw("colz SAME");
  //  pid_he4 ->Draw("colz SAME");
  //  c_dist -> Draw();
  //  m0->Draw();
  c1->SaveAs("padplane.png");
*/
  TCanvas *c2 = new TCanvas("c2","c2",1);
  //  pid_p->Draw("colz");
  c2->SaveAs("pid.png");
    
return 0;

}

