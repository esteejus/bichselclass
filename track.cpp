//track.cpp
#include "track.h"

void Track::SetMomentum(double mom){
  t_momentum=mom;
  bgamma=(mom/sqrt(mom*mom + t_mass*t_mass))*(sqrt(mom*mom + t_mass*t_mass)/t_mass);
  m_path = GetM0(bgamma);
  
    return;
}

void Track::Getfarray(){
  //probably should check if t_length is integer of x_seg
  double t_dist = 0.;//total distance traveled
  double eloss = 0.;
  double seg_dist = 0.;//current segment distance traveled
  double dx = 0.;//step size
  
  while(t_dist<=t_length){
    //      std::cout << "Current seg dist is "<<seg_dist<< " eloss is "<<eloss<<" current total l "<<t_dist<<std::endl;
    dx = GetMCstep();
    seg_dist += dx;
    t_dist += dx;
    
    if(seg_dist>=x_seg){
      //here we check if the next step dx, above, put us over our desired analyzed segment size
      seg_dist = seg_dist-x_seg;//remainder dx in next segment analyzed
      f_array.push_back(eloss/x_seg);
      eloss = GetEloss();//get the  energyloss in the next segment analyzed with current dx
      //	std::cout<<" dist is greater than seg "<<std::endl;
    }
      
    else eloss += GetEloss();//step size did not put us over the current analyzed segment; step in dE
      
  }


  return ;
}


TH2D* Track::GraphMomRange(int mc_tracks,int steps,double mom_min, double mom_max){
  double mom_step = (mom_max - mom_min)/steps;
  TString histname = Form("tracklength %f cm P10 gas",t_length);
  TH2D *pid = new TH2D("pid",histname,1000,0,5000,1000,0,20);
  pid->GetYaxis()->SetTitle("<C> [keV/cm]");
  pid->GetYaxis()->CenterTitle();
  pid->GetXaxis()->SetTitle("p [MeV/c]");
  pid->GetXaxis()->CenterTitle();
  
  for(int j=0;j<=steps;++j){
    for(int i=0;i<mc_tracks;++i){
      SetMomentum(mom_min+mom_step*j);
      Getfarray();
      SortArray();
      Truncate();
      GetC();
      for(int k=0;k<c_array.size();++k) pid->Fill(t_momentum,c_array.at(k));
      f_array.clear();
      c_array.clear();
    }
  }
 
 return pid;
}



TH1D* Track::Drawfdist(int mc_events,double max_eloss){
  TString histname = Form("test");
  TH1D *dist = new TH1D("f_dist",histname,1000,0,max_eloss);
  for(int i=0;i<mc_events;++i){
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

//void Track::SimpsonNInt(){

//}

//I want to get the array of f(xi) and xi array of the f(C) distribution
//this means i need to give the step size in xi and bin the f(C) values
//or I can just work from the TH1D but I think I dont want to do that.
//Once i have the two arrays then I can integrate the funciton
//I will make a function maybe SimponNInt() that an returns an array which is the integrated values
//of the funciton and the new x values of the integrand
//Then Within a program I can 

void Track::Truncate(){
  //factor is fraction  you keep. i.e. factor=.7 means we keep bottom 70% throw away top 30%
  int elem2trunc = round((1-t_factor)*f_array.size());
  //  std::cout << " elem to trunc is "<<elem2trunc<<std::endl;
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
