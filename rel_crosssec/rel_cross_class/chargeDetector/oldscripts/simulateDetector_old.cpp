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
#include "TCanvas.h"
#include "../dielectric.h"
#include "TRandom3.h"

using namespace std;

struct posvec {
  vector<double> x;
  vector<double> y;
  vector<double> z;
};

struct trackinfo{
  vector<double> init_pos;
  int num_seg;
  double seg_length;
  TGraph *dist;
};
  
void GenerateTrack(double z,double mass, double mom, double t_length, vector<vector<double>> &padplane);

Dielectric InitP10Gas()
{
 Dielectric argon{1.66e-3,38.,18,"argon"};//Argon @ 1 atm
  argon.SetPhotoCross("/home/justin/bichsel/bichselclass/rel_crosssec/argon_10_500_Sakamoto_ext.dat");

  //Set real and imaginary dielectric tables for faster computaiton
  if(argon.SetRealTable("/home/justin/bichsel/bichselclass/rel_crosssec/rel_cross_class/real_argon.dat") == false)
    {
      argon.GetRealDielectric();
      argon.WriteToFile("real");
    }
  if(argon.SetImgTable("/home/justin/bichsel/bichselclass/rel_crosssec/rel_cross_class/img_argon.dat") == false)
    {
      argon.GetImgDielectric();
      argon.WriteToFile("img");
    }  

  Dielectric ch4{.667e-3,16.,10,"ch4"};//CH4 @ 1 atm
  ch4.SetPhotoCross("/home/justin/bichsel/bichselclass/rel_crosssec/ch4.dat");

  //Set real and imaginary dielectric tables for faster computaiton
  if(ch4.SetRealTable("/home/justin/bichsel/bichselclass/rel_crosssec/rel_cross_class/real_ch4.dat") == false)
    {
      ch4.GetRealDielectric();
      ch4.WriteToFile("real");
    }
  if(ch4.SetImgTable("/home/justin/bichsel/bichselclass/rel_crosssec/rel_cross_class/img_ch4.dat") == false)
    {
      ch4.GetImgDielectric();
      ch4.WriteToFile("img");
    }  

Dielectric p10 = p10.MixGas(argon,ch4,.9,.1);
 p10.SetName("p10");
 cout<<"Density of "<<p10.GetName()<<" "<<p10.GetDensity()<<" Mol mass "<<p10.GetMolMass()<<" Z "<<p10.GetZ()<<endl;
  
  //Set real and imaginary dielectric tables for faster computaiton
  if(p10.SetRealTable("/home/justin/bichsel/bichselclass/rel_crosssec/rel_cross_class/real_p10.dat") == false)
    {
      p10.GetRealDielectric();
      p10.WriteToFile("real");
    }
  if(p10.SetImgTable("/home/justin/bichsel/bichselclass/rel_crosssec/rel_cross_class/img_p10.dat") == false)
    {
      p10.GetImgDielectric();
      p10.WriteToFile("img");
    }  


  return p10;
}

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

double GetFract(double x1, double x2, double y1, double y2, double x_0, double y_0)
{
  //x1 x2 integral in row direciton from x1 to x2 in current pad reference frame
  //y1 y2 integral in layer direciton from y1 to y2 in current pad reference frame
  double sigma_x = 3.4;
  double sigma_y = 3.4;
  //  double sigma_x = 5;
  //  double sigma_y = 5;

  double u1 = (x1 - x_0)/sigma_x;
  double u2 = (x2 - x_0)/sigma_x;
  double v1 = (y1 - y_0)/sigma_y;
  double v2 = (y2 - y_0)/sigma_y;

  return (1./4)*(erf(v2)-erf(v1))*(erf(u2)-erf(u1));
}

void FillPadPlane(double xpos, double ypos, double charge, vector<vector<double>> &padplane)
{


  int row = floor(xpos/8.);//current row
  int layer = floor(ypos/12.);//current layer

  double x_inpad = (xpos - row * 8.) - 4;//in coordiantes of tpc current pad
  double y_inpad = (ypos - layer * 12.) - 6;//in coordinates of tpc current pad

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


  double x_step = 8;
  double y_step = 12;
  int row_loop = 0;
  int layer_loop = 0;
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

	  
	  if( ((row + i) >= 0 && (layer + j) >= 0) && ((layer + j <= 112) && (row + i <= 108)) )
	    {
	      padplane.at(layer + j).at(row + i) += charge;//*GetFract(x1,x2,y1,y2,x_inpad,y_inpad);
	      //	      cout<<"curr row layer "<<row+i<<" "<<layer+j<<" "<< padplane.at(row + i).at(layer + j)<<endl;
	    }
	}
    }
  
  return;
}

/*
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
	  //	  cout<<"Filling charge "<<endl;
	  /// 	  cout<<"Row layer "<<row<<" "<<layer<<endl;
	  //	  cout<<"pad plane charge is "<<      padplane.at(layer).at(row)<<endl;
	  //	  cout<<"charge is "<<charge<<endl;
	}
      if(layer+1 <113 && row+1<109)
	padplane.at(layer + 1).at(row + 1) += charge;
    }
  
  return;
}
*/
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

double GetC(vector<vector<double>> padplane, double threshold, double trunc_f)
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
	  chg_layers.push_back(sum/1.2);//TOOK OUT /1.2 need if you are looking at energy loss
      //      cout<<"size "<<curLayer.size()<<endl;
    }

    std::sort(chg_layers.begin(),chg_layers.end(),[](double i, double j){return i < j;});
  int element = trunc_f*chg_layers.size();
  if(element != 0)
    element -= 1;//erase saves the elment you put in. to erase element a this pos sub 1
  chg_layers.erase(chg_layers.begin()+element,chg_layers.end());
  double sum = std::accumulate(chg_layers.begin(),chg_layers.end(),0.0);
  sum = sum/chg_layers.size();

  if(chg_layers.size()==0)
    sum = 0.;
  
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
	    if(iLayer == 2)cout<<"Plotting charge "<<iRow<<" charge "<<charge<<endl;
	    padplane_hist->SetBinContent(iLayer,iRow,charge);
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
  /*  double b_gamma = mom/mass;
  Bichsel b(b_gamma);
  b.SetInvXSec("P10M0invw_31623.inv");
  b.SetM0Table("m0.dat");

  double x_0 = 436;//initial x pos approx center of pad plane
  double y_0 = 0; //initial y pos really z axis on real TPC 
  double dist = 0; //distance traveled 
  while(dist < t_length)
    {
      double dx = b.GetMCstep()*10;//GetMCStep is in [cm]
      double dE = b.GetEloss();
      dist +=  dx;
      auto hit = GetCurrentPad(x_0,dist);//find hit from x,y info
      cout<<"Charge of event is "<<dE*pow(z,2)<<" "<<dist<<endl;
      auto hit_response = GetPadResponse(hit,dE*pow(z,2));
      FillPadPlane(hit_response,padplane);
    }
  */
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


/*
void GenerateTrackLinear(posvec &electrons, trackinfo t_info, TVector3 mom)
{
  electrons.x.clear();
  electrons.y.clear();
  electrons.z.clear();

  double x_curr = t_info.init_pos.at(0);
  double y_curr = t_info.init_pos.at(1);
  double z_curr = t_info.init_pos.at(2);
  int ns = t_info.num_seg;
  double seg_length = t_info.seg_length;
  TGraph *dist = t_info.dist;

  auto fcn = [=](Double_t *x, Double_t *) {return dist->Eval(x[0]);};
  double xmin, xmax, ymin,ymax;
  dist->GetPoint(0,xmin,ymin);
  dist->GetPoint(dist->GetN()-1,xmax,ymax);//index starts from 0
 
  TF1 *bichsel = new TF1("bichsel",fcn,xmin,xmax,0);//0 for number of parameters 

  double phi = atan2(mom.Z(),mom.X());
  double psi = atan2(mom.Y(),mom.Z());
  
  double stepsize = seg_length/(num_elec + 1); // add 1 for number of gaps between segment

  for(int iSeg = 0 ;iSeg < ns; iSeg++)
    {
      for( int iEl = 0 ;iEl < num_elec; iEl++)
	{
	  double pos_inseg = stepsize * iEl; 

	}
    }



}


	  double dr = ran -> Gaus(0,sigmaT);
	  double angle = ran -> Uniform(2*TMath::Pi());
	  double dx = dr * cos(angle);
	  double dz = dr * sin(angle);
	  

	  auto hit = GetCurrentPad(x_0,1);//find hit from x,z info
	  cout<<"Charge of event is "<<e_loss<<" "<<endl;
	  auto hit_response = GetPadResponse(hit,e_loss);
	  FillPadPlane(hit_response,padplane);

*/

int main()
{
  Dielectric p10 = InitP10Gas();

  int npoints = 1e3;
  double amu     = 931.5;     //MeV/c^2 from Nuclear Wallet cards
  int    ns      = 108; // number of segments
  double segment = 1.2; // [cm] segment analyzed
  double factor  = .7; // truncation factor

  //masses from nuclear wallet cards
  //masses from nuclear wallet cards
  double p_mass = amu+7.2889;
  double d_mass = amu*2+13.136;
  double t_mass = amu*3+14.95;
  double he3_mass = amu*3+14.931;
  double he4_mass = amu*4+2.425;
  
  //  double bgamma = .4;
  double p_mom = 903;
  double d_mom = 898;
  double t_mom = 1700;
  double he3_mom = 1700*2;
  double he4_mom = 1700*2;
    
  double part_mass = d_mass;
  double part_mom  = d_mom;


  double eVperelec = 26.4; //eV/elec

  vector<double> cdf_x,cdf_y, dist_x,dist_y;
  
  TGraph *  bichsel_p10 = p10.DrawBichselSeg(part_mass,part_mom,segment,npoints,0,14000e3);
  TH1D *  c_dist_p10 = p10.GetMCdist(bichsel_p10,"c_dist_p10",segment,ns,factor,1e3); 
  auto interp = [=](Double_t *x, Double_t *) {return bichsel_p10->Eval(x[0]);};
  TF1 *bichsel_fcn = new TF1("bichsel_fcn",interp,600,1e4,0);


  for(int iPoint = 0 ;iPoint < bichsel_p10->GetN(); iPoint++)
    {
      double x,y;
      bichsel_p10->GetPoint(iPoint,x,y);
      if(y > 1e-30)
	{
	  dist_x.push_back(x);
	  dist_y.push_back(y);
	} 
    }
  p10.GetCDF(dist_x,dist_y,cdf_x,cdf_y);

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *cdf_table = gsl_spline_alloc(gsl_interp_akima,cdf_x.size());
  gsl_spline_init(cdf_table,cdf_y.data(),cdf_x.data(),cdf_x.size());

  cout<<"integral expected "<<  c_dist_p10->Integral("width")<<endl;
  cout<<"mean expected "<<  c_dist_p10->GetMean()<<endl;
  /*
  TRandom3 *ran = new TRandom3(1234567);
  double garfield_coeff = .224;// mm/sqrt(cm);
  double diff_coeff = .071; // mm/sqrt(mm) transverse diffusion
  double drift_l = 200;// drift length [mm]
  double sigmaT = diff_coeff * sqrt(drift_l);

  TH1D * c_dist = new TH1D("c_dist","c_dist",1e3,0,1e4);


  cout<<"Starting MC section ==============="<<endl;
  int num_mc = 1e4;
  //  vector<vector<double>> padplane[num_mc];

    for(int iTrack = 0; iTrack < num_mc; iTrack++)
      {
	if(iTrack%1000==0)cout<<"Event "<<iTrack<<endl;
	double z_curr = 0; // init position
	double x_curr = 460; //int position

	vector<vector<double>> padplane = InitPadPlane();
	//  padplane = InitPadPlane();
      for(int iSeg = 0 ;iSeg < ns; iSeg++)
	{
	  double x = ran -> Uniform();
	  double e_loss =  gsl_spline_eval(cdf_table,x,acc);

	  int num_elec  = floor(e_loss/eVperelec);
	  double e_chg = (e_loss/eVperelec)/num_elec;

	  for( int iEl = 0 ;iEl < num_elec; iEl++)
	    {
	      z_curr = ran -> Uniform(iSeg *(segment *10), (iSeg + 1)*(segment*10));
	      //	      double dr = ran -> Gaus(0,sigmaT);
	      //	      double angle = x * 2 * TMath::Pi();//ran -> Uniform(2*TMath::Pi());
	      //	      double drift_x = dr * cos(angle) + x_curr;
	      //	      double drift_z = dr * sin(angle) + z_curr;
	      double drift_x =  x_curr;
	      double drift_z =  z_curr;
	      FillPadPlane(drift_x,drift_z,1.,padplane);
	    }
	}
      double c_eloss = GetC(padplane,.01,.7);
      c_dist->Fill(c_eloss*eVperelec-84);      
      }
  

 ofstream cfile;
  cfile.open("cstrag_theory_prf_diff.dat");
  
  for(int i = 1; i <= c_dist->GetXaxis()->GetNbins();i++)
    {
      double mean = c_dist->GetBinCenter(i);
      double value = c_dist->GetBinContent(i);

      if(value > 1e-10 && value > 0)
	cfile<<mean/1000<<"\t"<<value<<endl;
    }
  

        TH2D *pp = PlotPadPlane(padplane[1]);
  
  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->SetLogz();
  pp->GetYaxis()->SetRangeUser(55,60);
  pp->GetXaxis()->SetRangeUser(0,5);
  pp->Draw("colz");
  c2->SaveAs("padplane.png");
  
  */
  TCanvas *c1 = new TCanvas("c1","c1",1);
    c1->SetLogx();
  // bichsel_p10->Draw();
  //  c_dist->Scale(1./c_dist->Integral("width"));
    c_dist->GetXaxis()->SetRangeUser(4e3,7e3);
  //  c_dist->Draw();
  c_dist_p10->SetLineColor(2);
  c_dist_p10->Draw("same");

  c1->SaveAs("cdist_compare.png");

return 0;

}

