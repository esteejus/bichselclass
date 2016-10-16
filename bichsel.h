//bichsel.h
#ifndef BICHSEL_H
#define BICHSEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "TRandom3.h"
#include "TH1.h"
#include "TString.h"
#include "spline.h"

//using namespace std;

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

  void SetInvXSec (const std::string&);
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

#endif