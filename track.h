//trach.h
#ifndef TRACK_H
#define TRACK_H

#include "bichsel.h"
#include "TH2.h"
#include <iostream>

class Track : public Bichsel{
  double t_length,x_seg,t_factor, t_mass, t_momentum, bgamma;
  std::vector<double>c_array;//array for truncated mean distribution
  std::vector<double>f_array;//array for f(E) distribution
  
public:
  Track(double mass, double mom, double length,double seg, double factor) : t_mass(mass), t_momentum(mom), t_length(length),x_seg(seg), t_factor(factor),Bichsel((mom/sqrt(mom*mom + mass*mass))*(sqrt(mom*mom + mass*mass)/mass)){}
  
  TH1D* Drawfdist(int,double);
  TH1D* DrawCdist(double);

  TH2D* GraphMomRange(int,int,double,double);//graph a momentum range

  void SetLength(double length){t_length = length;}
  void SetTruncFactor(double factor){t_factor = factor;}
  void SortArray(){sort(f_array.begin(), f_array.end(), std::less<double>());}
  void Printfarray(){for(int i=0;i<f_array.size();++i) std::cout << f_array.at(i) << std::endl;}
  void PrintCarray(){for(int i=0;i<c_array.size();++i) std::cout << c_array.at(i) << std::endl;}
  void SetMomentum(double);
  void Getfarray();//f(E) is the probability distribution for an energy loss E; sample from this and store in f_array 
  void Truncate();
  void GetC();//calculate C for a given track
  void GetCArray(int);//get c for mc steps
  
  double GetMomentum(){return t_momentum;}
  double GetLength(){return t_length;}
  double GetTruncFactor(){return t_factor;}
  double GetCavg();
  double GetCsigma();
  double GetMean(TH1D *&dist){return dist->GetMean();}
  double GetFWHM(TH1D*&);
  double GetSigma(TH1D*&dist){return GetFWHM(dist)/2.355;}

  std::vector<std::vector<double>> HistArray(double,double,double);
  std::vector<std::vector<double>> SimpsonNInt(double,double,double);


};

#endif
