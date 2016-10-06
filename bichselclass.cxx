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

using namespace std;

bool wayToSort(double i, double j) { return i < j; }

class Bichsel {
protected:
  double m_path;//mean free path collisions/cm
  TRandom3 ran;
  std::vector<double>elosstable;
  std::vector<double>elossarray; //array to store energy loss simulations
public:
  Bichsel(double path) : m_path(path), ran(0) {}
  void SetInvXSec (const string&);
  void PrintInvTable();
  void GetElossArray(int);
  void SetArrayElem(double value){elossarray.push_back(value);}
  void PrintElossArray();
  void ClearArray(){elossarray.clear();}
  void SortArray(){sort(elossarray.begin(), elossarray.end(), wayToSort);}
  void Truncate(double);
  TH1D* DrawElossDist(int);
  TH1D* DrawMultColl(int,double);

  double GetMean();
  double GetEloss();
  double InterTableLin(double);
  double GetMCstep();
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
}

void Bichsel::PrintInvTable(){
  for(int i=0;i<elosstable.size();i++)cout<<i<<" "<<elosstable.at(i)<<endl;
  return;
}

double Bichsel::GetEloss(){
  double random=0;
  double r_eloss=-99;
  random = ran.Uniform(elosstable.size()-1);
  r_eloss = InterTableLin(random);
  
  return r_eloss;
}

//
void Bichsel::GetElossArray(int mc_events){
  double ran_array[mc_events];
  double random=0;
  ran.RndmArray(mc_events,ran_array);
  for(int i=0;i<mc_events;i++){
    //       random = ran_array[i]*(elosstable.size()-1);//uniform dist between elosstable range
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
  cout << " size of array "<<elossarray.size() <<endl;
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


class Track : public Bichsel{
  double t_length,x_seg;
  std::vector<double>c_array;//array for truncated mean
public:
  Track(double path,double length,double seg) :  Bichsel(path),t_length(length),x_seg(seg){}
  TH1D* GetfDist(int,double);//f(E) is the probability distribution for an energy loss E 
  //  void GetCdist();
  void SetMFree(double path){m_path=path;}
  void SetLength(double length){t_length=length;}
  double GetMFree(){return m_path;}
  double GetLength(){return t_length;}
  
};

TH1D* Track::GetfDist(int mc_events,double max_eloss){
  TH1D *dist = new TH1D("eloss","Monte Carlo Energy loss distribution",1000,0,max_eloss);
  //probably should check if t_length is integer of x_seg
  for(int i=0;i<mc_events;++i){
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
	//	SetArrayElem(eloss); //setting the eloss of the previous segment
	dist->Fill(eloss);//fill histigram with eloss of previous segment
	eloss = GetEloss();//get the  energyloss in the next segment analyzed with current dx
	//	cout<<" dist is greater than seg "<<endl;
      }
      
      else eloss += GetEloss();//step size did not put us over the current analyzed segment; step in dE
      
    }

    //    cout << "end of this track ------------------"<<endl<<endl;
    //    PrintElossArray();
    /*    SortArray();
	  Truncate(.7);
	  PrintElossArray();
    */
    //        ClearArray();

  }
  //dist->Scale(1.
  return dist;
}




int main(){
  //  Track t{40,10,2};
  //t.SetInvXSec("P10M0invw_31623.inv");
  //       t.GetElossArray(1e6);

  Bichsel b(40.);
  b.SetInvXSec("P10M0invw_31623.inv");
  b.GetElossArray(1e6);

  TH1D *h1;
  // h1 = b.DrawElossDist(100);
  //  h1 = t.GetfDist(1e6,20000);
  h1 = b.DrawMultColl(1,100);

  TCanvas *c1 = new TCanvas(1);
  h1->GetXaxis()->SetRangeUser(9,100);
  h1->Draw();
  c1->SetLogx();
  c1->SaveAs("test.jpg");

  return 0;
}
