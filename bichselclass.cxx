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
  double m_path;//mean free path collisions/cm
  TRandom3 ran;
  std::vector<double> elosstable;
  std::vector<double> elossarray; //array to store energy loss simulations
  std::vector<double> fvpx={0.316,0.398,0.501,0.631,0.794,1,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.31,7.943,10,12.589,15.849,19.953,25.119,31.623,39.811,50.119,63.096,79.433,100,125.893,158.489,199.526,251.189,316.228,398.107,501.187,630.957,794.328,1000};
  std::vector<double> fvpy={211.0726,146.5664,103.9873,76.0672,57.9161,46.2566,38.8999,34.3884,31.7545,30.357,29.7722,29.7206,30.018,30.543,31.2156,31.9825,32.8078,33.6658,34.5369,35.4067,36.2903,37.2469,38.055,38.6576,39.0968,39.4162,39.6515,39.8283,39.9648,40.0725,40.159,40.2288,40.285,40.3296,40.3634,40.3885}; 
  tk::spline fvp;
  //  fvp.set_points(fvpx,fvpy);  
public:
  Bichsel(double path) : m_path(path), ran(0) {}

  void SetInvXSec (const string&);
  void PrintInvTable();
  void GetElossArray(int);
  void SetArrayElem(double value){elossarray.push_back(value);}
  void PrintElossArray();
  void ClearArray(){elossarray.clear();}
  //  void SortArray(){sort(elossarray.begin(), elossarray.end(), wayToSort);}
  void Truncate(double);

  TH1D* DrawElossDist(int);
  TH1D* DrawMultColl(int,double);

  double GetMean();
  double GetEloss();
  double InterTableLin(double);
  double GetMCstep();
  double GetM0(double);
};

double Bichsel::InterTableLin(double ran){
  double y = -999999;//interpolated y value from table
  int rd_up=floor(ran);//round up value
  int rd_down=ceil(ran);//round down and down integers for 

  if(rd_up!=rd_down){
    double m=0.;//slope
    m=(elosstable.at(rd_up)-elosstable.at(rd_down))/(rd_up-rd_down);
    y=m*(ran-rd_down)+elosstable.at(rd_down);
  }
  //elosstable starts from 0 therefore size=rd_up (max) +1
  //error for array out of bounds
  else if(rd_down<0 || (rd_up+1)>elosstable.size()) cout << "ERROR : Interpolated point out of range of table" <<endl;

  else y=elosstable.at(rd_up); //ran must be an integer in double format
  
  return y;
}

void Bichsel::SetInvXSec (const string& filename) {
  ifstream file (filename.c_str());
  double d_value=0;
  string rnd_num;
  string index,line;

  if (file.is_open()){
    elosstable.push_back(0);// first value starts at .0001 not 0
    while(getline (file,line) ){
      
      istringstream in(line);
      in>>rnd_num;
      in>>d_value;
      elosstable.push_back(d_value);
    }
  }
  else cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << endl;

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
  fvp.set_points(fvpx,fvpy);

  return fvp(x);
}

double Bichsel::GetEloss(){
  double random=0;
  double r_eloss=-99;
  random = ran.Uniform(elosstable.size()-1);
  r_eloss = InterTableLin(random);
  
  return r_eloss;
}


void Bichsel::GetElossArray(int mc_events){
  double ran_array[mc_events];
  double random=0;
  ran.RndmArray(mc_events,ran_array);
  for(int i=0;i<mc_events;i++){
    random = ran.Uniform(elosstable.size()-1);//uniform dist between elosstable range
    elossarray.push_back(InterTableLin(random));
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
  double t_length,x_seg,t_factor;
  std::vector<double>c_array;//array for truncated mean distribution
  std::vector<double>f_array;//array for f(E) distribution
  
public:
  Track(double path,double length,double seg, double factor) :  Bichsel(path),t_length(length),x_seg(seg), t_factor(factor){}

  void Getfarray();//f(E) is the probability distribution for an energy loss E; sample from this and store in f_array 
  void SortArray(){sort(f_array.begin(), f_array.end(), wayToSort);}
  void Printfarray(){for(int i=0;i<f_array.size();++i) cout << f_array.at(i) << endl;}
  TH1D* Drawfdist(int,double);
  void GetCdist();
  void SetMFree(double path){m_path = path;}
  void SetLength(double length){t_length = length;}
  void SetTruncFactor(double factor){t_factor = factor;}
  void Truncate();

  //double GetC()DONT FORGET C = dE/dx !!!! check formula//return 
  double GetMFree(){return m_path;}
  double GetLength(){return t_length;}
  double GetTruncFactor(){return t_factor;}
  double GetMean(TH1D *&dist){return dist->GetMean();}
  double GetFWHM(TH1D*&);
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
      f_array.push_back(eloss);
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
    Getfarray();
    //        Printfarray();
	//     cout<<"after sort"<<endl;
    SortArray();
    // Printfarray();
    Truncate();
    //  Printfarray();
    //units of energy loss are in eV; which add up for multiple collisions in a segment
    //for 1cm segments the better unit is keV
for(int j=0;j<f_array.size();++j) dist->Fill(f_array.at(j)/1000);
    f_array.clear();
    //       Printfarray();
       //       cout<<"end =============="<<endl;
    //    cout<<"f array size"<<f_array.size();
  }
  dist->Scale(1./dist->GetEntries());

  return dist;
}

void Track::Truncate(){
  //factor is fraction  you keep. i.e. factor=.7 means we keep bottom 70% throw away top 30%
  int elem2trunc = round((1-t_factor)*f_array.size());
  //  cout << " elem to trunc is "<<elem2trunc<<endl;
  for(int i=0;i<elem2trunc;++i) f_array.pop_back();

  return;
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
  /*  std::vector<double> fvpx={0.316,0.398,0.501,0.631,0.794,1,1.259,1.585,1.995,2.512,3.162,3.981,5.012,6.31,7.943,10,12.589,15.849,19.953,25.119,31.623,39.811,50.119,63.096,79.433,100,125.893,158.489,199.526,251.189,316.228,398.107,501.187,630.957,794.328,1000};
    std::vector<double> fvpy={211.0726,146.5664,103.9873,76.0672,57.9161,46.2566,38.8999,34.3884,31.7545,30.357,29.7722,29.7206,30.018,30.543,31.2156,31.9825,32.8078,33.6658,34.5369,35.4067,36.2903,37.2469,38.055,38.6576,39.0968,39.4162,39.6515,39.8283,39.9648,40.0725,40.159,40.2288,40.285,40.3296,40.3634,40.3885};
    tk::spline fvp;

  fvp.set_points(fvpx,fvpy);
  cout<<fvp(.38)<<endl;
  */
  Track t{40,10,2,.7};
  t.SetInvXSec("P10M0invw_31623.inv");
  
  //       t.GetElossArray(1e6);
  //  Bichsel b(40.);
  //  b.SetInvXSec("P10M0invw_31623.inv");
  //  b.GetElossArray(1e6);
  //  std::vector<double> *x,*y;
   int steps =1e5;
    TGraph h1(steps);
    double scale=(5000-.2)/steps;
  for(int i=0;i<steps;++i){
    double xi = .2+scale*i;
    h1.SetPoint(i,xi,t.GetM0(xi));
  }
  /*    TH1D *h1;
    // h1 = b.DrawElossDist(100);
    //    h1 = t.Drawfdist(1e5,20);
    //    cout<<t.GetFWHM(h1)<<endl;
//  h1 = b.DrawMultColl(2,100);
*/
  TCanvas *c1 = new TCanvas(1);
  //  h1->GetXaxis()->SetRangeUser(9,100);

  h1.Draw();
    c1->SetLogx();
  c1->SaveAs("test.jpg");

  return 0;
}
