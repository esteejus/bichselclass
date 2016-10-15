//bichsel.cpp
#include "bichsel.h"

double Bichsel::InterTableLin(double r){
  double y = -999999;//interpolated y value from table
  int rd_down=floor(r);//round up value
  int rd_up=ceil(r);//round down and down integers for 

  if(rd_up!=rd_down){
    double m=0.;//slope
    m=(elosstable.at(rd_up)-elosstable.at(rd_down))/(rd_up-rd_down);
    y=m*(r-rd_down)+elosstable.at(rd_down);
  }
  //elosstable starts from 0 therefore size=rd_up (max) +1
  //error for array out of bounds
  else if(rd_down<0 || (rd_up+1)>elosstable.size()) std::cout << "ERROR : Interpolated point out of range of table" <<std::endl;

  else y=elosstable.at(rd_up); //ran must be an integer in double format
  
  return y;
}

double Bichsel::InterTableSpline(double r){
  if(r<0){
    std::cout << "ERROR rand num is less than 0";
    std::cout << "Please check random number code limits" << std::endl;
    //DO I NEED TO EXIT!!!???
  }
  else if(r>1){
    std::cout << "Random number is " << r <<std::endl;
    std::cout << "ERROR rand is greater than 1 which is max table value" << std::endl;
    std::cout << "Please check random number code limits" << std::endl;
    //EXIT HERE TOO!!!!
  }

  return table(r);
}

void Bichsel::SetInvXSec (const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_ran=0;
  //  string rnd_num;
  std::string index,line;

  if (file.is_open()){
    //    elosstable.push_back(0);// first value starts at .0001 not 0
    while(getline (file,line) ){
      
      std::istringstream in(line);
      in>>d_ran;
      in>>d_value;
      elosstable.push_back(d_value);
      tablex.push_back(d_ran);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;
  table.set_points(tablex,elosstable);

  return;
}

void Bichsel::PrintInvTable(){
  for(int i=0;i<elosstable.size();i++)std::cout<<i<<" "<<elosstable.at(i)<<std::endl;

  return;
}

double Bichsel::GetM0(double x){
  if(x<.25){
    std::cout << "Beta*gamma is less than the min. range for interpolation" << std::endl;
    std::cout << "setting Beta*gamma to a reasonable minimum .25" << std::endl;
    x=.25;
  }
  else if(x>5000){
    std::cout << "Beta*gamma is less than the max range for interpolation" << std::endl;
    std::cout << "setting Beta*gamma to a reasonable maximum 5000" << std::endl;
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
  std::cout << "Energy loss array is in [keV] "<<std::endl;
  for(int i=0;i<elossarray.size();++i)std::cout << elossarray.at(i) <<std::endl;

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
