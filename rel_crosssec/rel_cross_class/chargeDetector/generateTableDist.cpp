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
#include "../dielectric.h"
#include "TRandom3.h"

using namespace std;

vector<double> erf_x,erf_y;
  double erf_min = -10;//set by min and max layer/ row ever seen
  double erf_max = 10;
  int num_p_erf = 1e4;

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

void InitERFTable()
{
  double step_size = (erf_max - erf_min)/num_p_erf;

  for(int i = 0; i < num_p_erf; i++)
    {
      double x = step_size * i + erf_min;
      double y = erf(x);
      erf_x.push_back(x);
      erf_y.push_back(y);
    }

  return;
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

  double u1 = (x1 - x_0)/(sqrt(2)*sigma_x);
  double u2 = (x2 - x_0)/(sqrt(2)*sigma_x);
  double v1 = (y1 - y_0)/(sqrt(2)*sigma_y);
  double v2 = (y2 - y_0)/(sqrt(2)*sigma_y);

  double step_size = (erf_max - erf_min)/num_p_erf;

  int idx = (v2 - erf_min)/step_size;
  double erf_v2 = erf_y.at(idx);

  idx = (v1 - erf_min)/step_size;
  double erf_v1 = erf_y.at(idx);

  idx = (u1 - erf_min)/step_size;
  double erf_u1 = erf_y.at(idx);

    idx = (u2 - erf_min)/step_size;
  double erf_u2 = erf_y.at(idx);

  //  return (1./4)*(erf(v2)-erf(v1))*(erf(u2)-erf(u1));
  return (1./4)*(erf_v2-erf_v1)*(erf_u2-erf_u1);
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
	      padplane.at(layer + j).at(row + i) += charge*GetFract(x1,x2,y1,y2,x_inpad,y_inpad);
	      //	      cout<<"curr row layer "<<row+i<<" "<<layer+j<<" "<< padplane.at(row + i).at(layer + j)<<endl;
	    }
	}
    }
  
  return;
}

double GetC(vector<vector<double>> padplane, double threshold, double trunc_f)
{
  int num_layers = 112;
  int num_rows = 108;
  vector<double> chg_layers = {};
  for(int iLayer = 1; iLayer<= num_layers; iLayer++)
    {
      padplane.at(iLayer).erase(remove_if(padplane.at(iLayer).begin(),padplane.at(iLayer).end(),[threshold](double x){ return x<threshold;}),padplane.at(iLayer).end());

      double sum = std::accumulate(padplane.at(iLayer).begin(),padplane.at(iLayer).end(),0.0);
      if(sum > 0)
	  chg_layers.push_back(sum/1.2);//TOOK OUT /1.2 need if you are looking at energy loss
      //      cout<<"size "<<padplane.at(iLayer).size()<<endl;
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


TH1D * FoldDetectorEffects(double part_mass, double part_mom, int ns, double segment, int npoints)
{
  Dielectric p10 = InitP10Gas();
  double eVperelec = 26.4; //eV/elec
  vector<double> cdf_x,cdf_y, dist_x,dist_y;

  TGraph *  bichsel_p10 = p10.DrawBichselSeg(part_mass,part_mom,segment,npoints,0,14000e3);
  TH1D * folded_dist = new TH1D("folded_dist","folded_dist",1e5,0,1e6);
  
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

  TRandom3 *ran = new TRandom3(1234567);
  double garfield_coeff = .224;// mm/sqrt(cm);
  double diff_coeff = .071; // mm/sqrt(mm) transverse diffusion
  double drift_l = 200;// drift length [mm]
  double sigmaT = diff_coeff * sqrt(drift_l);

  cout<<"Starting MC section ==============="<<endl;
  int num_mc = 1e5;
  //  vector<vector<double>> padplane[num_mc];

    for(int iTrack = 0; iTrack < num_mc; iTrack++)
      {
	auto padplane = InitPadPlane();
	if(iTrack%1000==0)cout<<"Event "<<iTrack<<endl;
	double z_curr = 0; // init position
	double x_curr = 460; //int position

	//	auto padplane = InitPadPlane();
      for(int iSeg = 0 ;iSeg < ns; iSeg++)
	{
	  double z_curr = iSeg * (segment* 10);
	  double x = ran -> Uniform();
	  double e_loss =  gsl_spline_eval(cdf_table,x,acc);

	  int num_elec  = floor(e_loss/eVperelec);
	  int rem_elec = num_elec%3;//remainder electrons
	  int split_chg = (num_elec - rem_elec)/3;

	  //for the three anode wires in a pad
	  FillPadPlane(x_curr, z_curr + 2 , split_chg, padplane);
	  FillPadPlane(x_curr, z_curr + 6 , split_chg, padplane);
	  FillPadPlane(x_curr, z_curr + 10, split_chg, padplane);

	  // cout<<"number of remainder electrons "<<rem_elec<<endl;
	  for(int iRem = 0; iRem < rem_elec; iRem++)
	    {
	      double gap = ran -> Uniform();
	      if(gap < (1./3))
		FillPadPlane(x_curr, z_curr + 2 , 1., padplane);
	      else if( gap < (2./3) && gap >= (1./3) )
		FillPadPlane(x_curr, z_curr + 2 , 1., padplane);
	      else if( gap < 1. && gap >= (2./3) )
		FillPadPlane(x_curr, z_curr + 2 , 1., padplane);
	    }
	  
	}

      double c_eloss = GetC(padplane,.01,.7);
      folded_dist->Fill(c_eloss*eVperelec);      

      }

    gsl_interp_accel_free(acc);
    gsl_spline_free(cdf_table);

    return folded_dist;
}

int main()
{

  InitERFTable();

  int npoints = 1e3;
  double amu     = 931.5;     //MeV/c^2 from Nuclear Wallet cards
  int    ns      = 108; // number of segments
  double segment = 1.2; // [cm] segment analyzed
  double factor  = .7; // truncation factor

  //masses from nuclear wallet cards
  //masses from nuclear wallet cards

  int nbg_points = 5e1;
  double bg_l = .1;
  double bg_h = 3;
  double bg_step = (log(bg_h)-log(bg_l))/nbg_points;
  double bg = 0;
  
  double part_mass = amu*2+13.136;
  double part_mom  = bg * part_mass;
  TH1D *hist = nullptr;
  TH1D *cdist[nbg_points];

  TFile *root_out = new TFile("inputtable.root","RECREATE");
  TTree *table = new TTree("table","table");
  table->Branch("bg",&bg);
  table->Branch("cdist",&hist);

  for(int iBg = 0; iBg < nbg_points; iBg++)
    {
      bg = bg_step * iBg + log(bg_l);
      bg = exp(bg);
      part_mom = bg * part_mass;
      cout<<"Bg is "<<bg<<endl;
      cdist[iBg] = FoldDetectorEffects(part_mass,part_mom,ns,segment,npoints);
      TString name = Form("bg_%f_ns_%i",bg,ns);
      cdist[iBg]->SetName(name);
      cdist[iBg]->Scale(1./cdist[iBg]->Integral("width"));
      hist = cdist[iBg];
      table->Fill();
    }
  table->Write();
  root_out->Close();
  
return 0;

}

