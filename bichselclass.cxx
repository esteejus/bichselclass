//class for bichsel straggling funcitons
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>
#include "TRandom3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TGraph.h"
#include "spline.h"

using namespace std;

bool wayToSort(double i, double j) { return i < j; }

class Bichsel {
protected:
  double bgamma,m_path;//mean free path collisions/cm
  TRandom3 ran;
  std::vector<double> tablex; //array to store x values (rand num 0-1) from Cumulative dist of eloss
  std::vector<double> elosstable; //array to store values of cumulative dist
  std::vector<double> elossarray; //array to store energy loss simulations

  //The arrays below come from Hans Bichsel's table ??? for P10 gas
  // fvpx represents the beta gamma values in this table for FVP theory
  // fvpy represents the M0 values [collisions/cm] for FVP theory
  std::vector<double> fvpx={0.316,0.398,0.501,0.631,0.794,1,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.31,7.943,10,12.589,15.849,19.953,25.119,31.623,39.811,50.119,63.096,79.433,100,125.893,158.489,199.526,251.189,316.228,398.107,501.187,630.957,794.328,1000};
  std::vector<double> fvpy={211.0726,146.5664,103.9873,76.0672,57.9161,46.2566,38.8999,34.3884,31.7545,30.357,29.7722,29.7206,30.018,30.543,31.2156,31.9825,32.8078,33.6658,34.5369,35.4067,36.2903,37.2469,38.055,38.6576,39.0968,39.4162,39.6515,39.8283,39.9648,40.0725,40.159,40.2288,40.285,40.3296,40.3634,40.3885}; 
  tk::spline fvp;
  tk::spline table;//eloss taable for spline
public:
  Bichsel(double bg) : bgamma(bg), ran(0){fvp.set_points(fvpx,fvpy); m_path = GetM0(bg);}

  void SetInvXSec (const string&);
  void PrintInvTable();
  void GetElossArray(int);
  void SetArrayElem(double value){elossarray.push_back(value);}
  void Setbg(double bg){bgamma = bg;}
  void PrintElossArray();
  void ClearArray(){elossarray.clear();}
  //  void SortArray(){sort(elossarray.begin(), elossarray.end(), wayToSort);}
  void Truncate(double);

  TH1D* DrawElossDist(int);
  TH1D* DrawMultColl(int,double);

  double GetMean();
  double GetEloss();
  double Getbg(){return bgamma;}
  double GetMpath(){return m_path;}
  double InterTableLin(double);
  double InterTableSpline(double);
  double GetMCstep();
  double GetM0(double);
};

double Bichsel::InterTableLin(double r){
  double y = -999999;//interpolated y value from table
  int rd_up=floor(r);//round up value
  int rd_down=ceil(r);//round down and down integers for 

  if(rd_up!=rd_down){
    double m=0.;//slope
    m=(elosstable.at(rd_up)-elosstable.at(rd_down))/(rd_up-rd_down);
    y=m*(r-rd_down)+elosstable.at(rd_down);
  }
  //elosstable starts from 0 therefore size=rd_up (max) +1
  //error for array out of bounds
  else if(rd_down<0 || (rd_up+1)>elosstable.size()) cout << "ERROR : Interpolated point out of range of table" <<endl;

  else y=elosstable.at(rd_up); //ran must be an integer in double format
  
  return y;
}

double Bichsel::InterTableSpline(double r){
  if(r<0){
    cout << "ERROR rand num is less than 0";
    cout << "Please check random number code limits" << endl;
    //DO I NEED TO EXIT!!!???
  }
  else if(r>1){
    cout << "Random number is " << r <<endl;
    cout << "ERROR rand is greater than 1 which is max table value" << endl;
    cout << "Please check random number code limits" << endl;
    //EXIT HERE TOO!!!!
  }

  return table(r);
}

void Bichsel::SetInvXSec (const string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_ran=0;
  //  string rnd_num;
  string index,line;

  if (file.is_open()){
    //    elosstable.push_back(0);// first value starts at .0001 not 0
    while(getline (file,line) ){
      
      istringstream in(line);
      in>>d_ran;
      in>>d_value;
      elosstable.push_back(d_value);
      tablex.push_back(d_ran);
    }
  }
  else cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << endl;
  table.set_points(tablex,elosstable);

  return;
}

void Bichsel::PrintInvTable(){
  for(int i=0;i<elosstable.size();i++)cout<<i<<" "<<elosstable.at(i)<<endl;

  return;
}

double Bichsel::GetM0(double x){
  if(x<.25){
    cout << "Beta*gamma is less than the min. range for interpolation" << endl;
    cout << "setting Beta*gamma to a reasonable minimum .25" << endl;
    x=.25;
  }
  else if(x>5000){
    cout << "Beta*gamma is less than the max range for interpolation" << endl;
    cout << "setting Beta*gamma to a reasonable maximum 5000" << endl;
    x=5000;
  }
  //    fvp.set_points(fvpx,fvpy);

  return fvp(x);
}

double Bichsel::GetEloss(){
  double random=0;
  double r_eloss=-99;
  //    random = ran.Uniform(elosstable.size()-1);
  //    r_eloss = InterTableLin(random);
  random = ran.Rndm();
  r_eloss = InterTableSpline(random);
  
  return r_eloss;
}


void Bichsel::GetElossArray(int mc_events){
  double ran_array[mc_events];
  double random=0;
  ran.RndmArray(mc_events,ran_array);
  for(int i=0;i<mc_events;i++){
    //    random = ran.Uniform(elosstable.size()-1);//uniform dist between elosstable range
    //    elossarray.push_back(InterTableLin(random));
    random = ran.Rndm();
    elossarray.push_back(InterTableSpline(random));
  }

  return;
}

void Bichsel::PrintElossArray(){
  cout << "Energy loss array is in [keV] "<<endl;
  for(int i=0;i<elossarray.size();++i)cout << elossarray.at(i) <<endl;

  return;
}

TH1D* Bichsel::DrawElossDist(int max_eloss=1e6){
  TH1D *dist = new TH1D("eloss","Monte Carlo Energy loss distribution",1000,0,max_eloss);
  for(int i=0;i<elossarray.size();++i) dist->Fill(elossarray.at(i));
  dist->Scale(1./elossarray.size());

  return dist;
}

TH1D* Bichsel::DrawMultColl(int num_coll,double max_eloss){

  TString filename = Form("Distribution for %i collisions",num_coll);
  TH1D *dist = new TH1D("eloss",filename,1000,0,max_eloss);
  int mc_events=1e6;
  double eloss=0;
  for(int i=0;i<mc_events;++i){
    eloss=0;
    for(int j=0;j<num_coll;++j) eloss += GetEloss();
    dist->Fill(eloss);
  }
  dist->Scale(1./(num_coll*mc_events));

  return dist;
}


double Bichsel::GetMCstep(){
  double dx=0;//dx step given by mean free path
  dx=-log(ran.Rndm())/m_path;

  return dx;
}

/*
  void Bichsel::Truncate(double factor){
  //factor is fraction  you keep. i.e. factor=.7 means we keep bottom 70% throw away top 30%
  int elem2trunc = round((1-factor)*elossarray.size());
  cout << " elem to trunc is "<<elem2trunc<<endl;
  for(int i=0;i<elem2trunc;++i) elossarray.pop_back();

  return;
  }

  double Bichsel::GetMean(){
  double sum=0.;
  for(int i=0;i<elossarray.size();++i) sum+=elossarray.at(i);

  return sum/elossarray.size();
  }
*/

class Track : public Bichsel{
  double t_length,x_seg,t_factor, t_mass, t_momentum, bgamma;
  std::vector<double>c_array;//array for truncated mean distribution
  std::vector<double>f_array;//array for f(E) distribution
  
public:
  //  Track(double path,double length,double seg, double factor) :  Bichsel(path),t_length(length),x_seg(seg), t_factor(factor){}
  Track(double mass, double mom, double length,double seg, double factor) : t_mass(mass), t_momentum(mom), t_length(length),x_seg(seg), t_factor(factor),Bichsel((mom/sqrt(mom*mom + mass*mass))*(sqrt(mom*mom + mass*mass)/mass)){}
  
  TH1D* Drawfdist(int,double);
  TH1D* DrawCdist(double);

  void Getfarray();//f(E) is the probability distribution for an energy loss E; sample from this and store in f_array 
  void SortArray(){sort(f_array.begin(), f_array.end(), wayToSort);}
  void Printfarray(){for(int i=0;i<f_array.size();++i) cout << f_array.at(i) << endl;}
  void SetMomentum(double mom){t_momentum=mom; Setbg((mom/sqrt(mom*mom + t_mass*t_mass))*(sqrt(mom*mom + t_mass*t_mass)/t_mass));}
  void SetLength(double length){t_length = length;}
  void SetTruncFactor(double factor){t_factor = factor;}
  void Truncate();
  void GetC();//calculate C for a given track

  double GetCavg();
  double GetCsigma();
  double GetLength(){return t_length;}
  double GetTruncFactor(){return t_factor;}
  double GetMean(TH1D *&dist){return dist->GetMean();}
  double GetFWHM(TH1D*&);
  double GetSigma(TH1D*&dist){return GetFWHM(dist)/2.355;}
};


void Track::Getfarray(){
  //probably should check if t_length is integer of x_seg
  double t_dist = 0.;//total distance traveled
  double eloss = 0.;
  double seg_dist = 0.;//current segment distance traveled
  double dx = 0.;//step size
  
  while(t_dist<=t_length){
    //      cout << "Current seg dist is "<<seg_dist<< " eloss is "<<eloss<<" current total l "<<t_dist<<endl;
    dx = GetMCstep();
    seg_dist += dx;
    t_dist += dx;
    
    if(seg_dist>=x_seg){
      //here we check if the next step dx, above, put us over our desired analyzed segment size
      seg_dist = seg_dist-x_seg;//remainder dx in next segment analyzed
      f_array.push_back(eloss/x_seg);
      eloss = GetEloss();//get the  energyloss in the next segment analyzed with current dx
      //	cout<<" dist is greater than seg "<<endl;
    }
      
    else eloss += GetEloss();//step size did not put us over the current analyzed segment; step in dE
      
  }


  return ;
}

TH1D* Track::Drawfdist(int mc_events,double max_eloss){
  TString histname = Form("test");
  TH1D *dist = new TH1D("f_dist",histname,1000,0,max_eloss);
  for(int i=0;i<mc_events;++i){
    if(i%1000==0)cout<<i<<endl;
    Getfarray();
    SortArray();
    Truncate();
    GetC();
    //units of energy loss are in eV/cm; which add up for multiple collisions in a segment
    //the better unit is keV/cm
    for(int j=0;j<f_array.size();++j) dist->Fill(f_array.at(j)/1000);
    f_array.clear();
  }
  dist->Scale(1./dist->GetEntries());

  return dist;
}

TH1D * Track::DrawCdist(double max_eloss){
  TString histname = Form("test");
  TH1D *dist = new TH1D("c_dist",histname,1000,0,max_eloss);
  for(int i = 0;i<c_array.size();++i) dist->Fill(c_array.at(i));

  return dist;
}

void Track::Truncate(){
  //factor is fraction  you keep. i.e. factor=.7 means we keep bottom 70% throw away top 30%
  int elem2trunc = round((1-t_factor)*f_array.size());
  //  cout << " elem to trunc is "<<elem2trunc<<endl;
  for(int i=0;i<elem2trunc;++i) f_array.pop_back();

  return;
}

void Track::GetC(){
  double sum=0.;
  for(int i=0;i<f_array.size();++i){
    sum += f_array.at(i)/1000;//scale to [keV/cm] from [eV/cm]
  }
  c_array.push_back(sum/f_array.size());

  return;
}

double Track::GetCavg(){
  double sum = 0.;
  for(int i=0;i<c_array.size();++i){
    sum += c_array.at(i);
  }
  
  return (sum/c_array.size());
}

double Track::GetCsigma(){
  double avg = GetCavg(); //<C> 
  double sigma =0.;
  for(int i=0;i<c_array.size();++i){
    sigma += (c_array.at(i) - avg)*(c_array.at(i) - avg)/c_array.size();
    //    cout << "c is "<< c_array.at(i)<<endl;
  }

  return sqrt(sigma);
}

double Track::GetFWHM(TH1D*&dist){
  int l_bin,u_bin;
  double mean = dist->GetMean();
  int mean_bin = dist->GetXaxis()->FindBin(mean);
  double value = dist->GetBinContent(mean_bin);//distribution value at mean position
  double fwhm = 0;
  l_bin = dist->FindFirstBinAbove(value/2);
  u_bin = dist->FindLastBinAbove(value/2);
  fwhm = dist->GetBinCenter(u_bin) - dist->GetBinCenter(l_bin);
  
  return fwhm;
}

//make a class called particle
//give the mass,momentum,
//calculate Betagamma, 
//later do for different z

int main(){
  double length = 40; // [cm] length of track
  double segment = 2; // [cm] segment analyzed
  double factor = .6; // truncation factor
  
  Track pion{140,56,length,segment,factor};

  pion.SetInvXSec("P10M0invw_31623.inv");
  cout << "Bgamma is " << pion.Getbg() << endl;
  cout << "Mean free path " << pion.GetMpath() <<endl;
  //  pion.SetMomentum(200);
  cout << "Bgamma is " << pion.Getbg() << endl;
  cout << "Mean free path " << pion.GetMpath() <<endl;
  //       t.GetElossArray(1e6);
  //What i need to do next
  /* 
I need to make sure the truncation works
the Draw fDist is not C dist
The C distribution is right in units
make a <C> function
verify the distribution f(C) works in han's paper
figure out how to make a graph of Sigma vs momentum for particle type

  */
  TH1D *h1;
  //   h1 = b.DrawElossDist(100);
  h1 = pion.Drawfdist(1e3,20);
  //  cout << "FWHM        " << pion.GetFWHM(h1) << endl;
  cout << "Mean: " << pion.GetCavg() <<" Sigma: " << pion.GetCsigma()<<endl;
  //        h1 = t.DrawMultColl(2,100);

  TCanvas *c1 = new TCanvas(1);
  //  h1->GetXaxis()->SetRangeUser(9,100);
  
  h1->Draw();
  //  c1->SetLogx();
  c1->SaveAs("test.jpg");

  return 0;
}
