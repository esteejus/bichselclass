//track.cpp
#include "track.h"

//things left to do
//Find a better way to call the distributions
//make sure that after calling the dist the arrays are cleared
//Check if you call two differenct instances of the class the TH1 is not rewrittedn for some reason it seem to be that way
//check that the track length is an integer of the segment size
//need to write a program to scale/fit the cumulative integrals

void Track::SetMomentum(double mom){
  t_momentum=mom;
  bgamma=(mom/sqrt(mom*mom + t_mass*t_mass))*(sqrt(mom*mom + t_mass*t_mass)/t_mass);
  m_path = InterpolateM0(bgamma);
  
    return;
}

void Track::Getfarray(){
    f_array.clear();
  //probably should check if t_length is integer of x_seg
  double t_dist = 0.;//total distance traveled
  double eloss = 0.;
  double seg_dist = 0.;//current segment distance traveled
  double dx = 0.;//step size

  while(t_dist<=t_length){
    dx = GetMCstep();
    seg_dist += dx;
    t_dist += dx;
    if(seg_dist>=x_seg){
      //here we check if the next step dx, above, put us over our desired analyzed segment size
      seg_dist = seg_dist-x_seg;//remainder dx in next segment analyzed
      f_array.push_back(eloss/x_seg);
      eloss = GetEloss();//get the  energyloss in the next segment analyzed with current dx
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
  TH1D *dist = new TH1D("f_dist",histname,2000,0,max_eloss);
  for(int i=0;i<mc_events;++i){
    Getfarray();
    SortArray();
    Truncate();
    GetC();
    //units of energy loss are in eV/cm; which add up for multiple collisions in a segment
    //the better unit is keV/cm
    int size_a = f_array.size();
    for(int j=0;j<size_a;++j) dist->Fill(f_array.at(j)/1000);
    f_array.clear();
  }
  dist->Scale(1./dist->Integral("width"));
  std::cout<<"Integral is "<<dist->Integral();
  return dist;
}

TH1D * Track::DrawCdist(double max_eloss){
  TString histname = Form("test");
  TH1D *dist = new TH1D("c_dist",histname,1000,0,max_eloss);
  int size_a = c_array.size();
  for(int i = 0;i<size_a;++i) dist->Fill(c_array.at(i));

  return dist;
}

std::vector<std::vector<double>> Track::HistArray(double step_size,double low_lim, double up_lim){

  int num_steps = ceil((up_lim-low_lim)/step_size);//round up value
  std::vector<std::vector<double> > hist;
  
  //sort the <C> array to prepare it for binning into vector
  sort(c_array.begin(),c_array.end(),std::less<double>());

  //Initialize an array integral[i][j] where i runs from 0,1
  //0 componet stores the  center of the bin value [keV/cm]
  //1 componet stores the counts of the histogram or normalized value
  hist.resize(2);

  for (int i = 0; i < 2; ++i) {
    hist[i].resize(num_steps);
  }
  
  //first we bin into a histogram
  for(int i=0;i<num_steps;++i){
    double center_x = (low_lim + step_size/2) + (step_size*i);
    int count=0;
    int size_a = c_array.size();
    for(int j=0;j<size_a;++j){
      if(c_array.at(j)>=(i*step_size) && ((i+1)*step_size)>=c_array.at(j)) count++;
    }
      hist[0][i] = center_x;
      hist[1][i] = count;
  }

  return hist;
}

std::vector<std::vector<double>> Track::SimpsonNInt(double step_size,double low_lim, double up_lim){

  int num_steps = ceil((up_lim-low_lim)/step_size);//round up value

  std::vector<std::vector<double> > integral,hist;
  hist = HistArray(step_size,low_lim,up_lim);

  //Initialize an array integral[i][j] where i runs from 0,1
  //0 componet stores the  x value in the cumulative integral [keV/cm]
  //1 componet stores the integral value which is normalized
  integral.resize(2);

  if(hist[0].size()!=num_steps)std::cout << "ERROR histogram and step size are different" << std::endl;

  if(hist[0].size()%2==0){std::cout<<"even"<<std::endl;hist[0].pop_back(); hist[1].pop_back();}//simpson's rule requies odd number set

  double sum = 0.;//cumulative integral
  double x = hist[0][0];//cumulative x posiiton
  
  for(int i=0;i<num_steps;i+=2){
    double dx =0.;//cumulative x position
    double df =0.;//integral step
    if( (i+2) != num_steps ){
      df = (step_size/3)*(hist[1][i] + 4*hist[1][i+1] + hist[1][i+2]);
      dx = (hist[0][i+2]-hist[0][i]);
      
      sum+=df;
      x+=dx;
      
      integral[0].push_back(x);
      integral[1].push_back(sum);
    }
    
    else break;
  }

  //last point is not a real physical point
  integral[0].pop_back();
  integral[1].pop_back();

  //scale the cumulative integral to 1
  for(int i=0;i<integral[0].size();++i) integral[1].at(i)=integral[1].at(i)/integral[1].back();

  return integral;
}

void Track::Truncate(){
  //factor is fraction  you keep. i.e. factor=.7 means we keep bottom 70% throw away top 30%
  int elem2trunc = round((1-t_factor)*f_array.size());
  for(int i=0;i<elem2trunc;++i) f_array.pop_back();

  return;
}

void Track::GetCArray(int nevents){
  c_array.clear();
  for(int i = 0;i < nevents;++i)
    {
      if(i%1000==0)std::cout<<"On event "<<i<<std::endl;
      f_array.clear();
      Getfarray();
      SortArray();
      Truncate();
      GetC();
    }
  return;
}

void Track::GetC(){
  double sum=0.;
  int size_a = f_array.size();
  for(int i=0;i<size_a;++i){
    sum += f_array.at(i)/1000;//scale to [keV/cm] from [eV/cm]
  }
  c_array.push_back(sum/f_array.size());

  return;
}

double Track::GetCavg(){
  double sum = 0.;
  int size_a = c_array.size();
  for(int i=0;i<size_a;++i){
    sum += c_array.at(i);
  }
  
  return (sum/c_array.size());
}

double Track::GetCsigma(){
  double avg = GetCavg(); //<C> 
  double sigma =0.;
  int size_a = c_array.size();
  for(int i=0;i<size_a;++i){
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
