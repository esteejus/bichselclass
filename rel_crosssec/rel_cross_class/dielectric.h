#ifndef DIELECTRIC_H
#define DIELECTRIC_H

#include <gsl/gsl_spline.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <functional>
#include <gsl/gsl_integration.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include <algorithm>

using namespace std;

class gsl_function_pp : public gsl_function
 {
    private:
   std::function<double(double)> _func;
    static double invoke(double x, void *params) {
     return static_cast<gsl_function_pp*>(params)->_func(x);
   }

    public:
 gsl_function_pp(std::function<double(double)> const& func) : _func(func){
      function=&gsl_function_pp::invoke;
      params=this;
    }    
};

class Dielectric {
  
 protected:
  double energy_p = 0;
  double avogadro    =  6.022e23;
  double m_elec      =  5.11e5;//mass of electron ev/c^2
  double hbar_c      =  1.97327e-5; // [eV *cm]
  double units_cross =  1e-18; //units the cross section is given in in cm^2
  double density = -1 , molarmass = -1 , atom_cm3 = -1;
  double z = -1; //z is atomic number
  string name;
  int moment = 0;
  
  std::vector<double> photoenergy; //array to store values of cumulative dist
  std::vector<double> photovalue; //array to store values of cumulative dist

  std::vector<double> img_energy, img_value;
  std::vector<double> real_energy, real_value;

  std::vector<double> relcross_energy, relcross_value;
  
  bool set_table = false;
  bool set_img   = false;//if img is calculated or set through table this is TRUE
  bool set_real  = false;//if real is calculated or set through table this is TRUE
  bool set_rel   = false;//if relatavistic cross section is set this is TRUE 

 public:
 Dielectric( double d_density, double d_molarmass, double zz, string nn) : density(d_density), molarmass(d_molarmass), z(zz), name(nn){ atom_cm3 = (d_density/d_molarmass) *avogadro;}
  
  string GetName(){return name;}

  void SetPhotoCross ( const std::string &);
  bool SetRealTable ( const std::string &);
  bool SetImgTable ( const std::string &);
  void SetDensity(double den){density = den;}
  void SetMolMass(double mm){molarmass = mm;}
  void SetAtomicNum(double zz){z = zz;}

  double GetDensity(){ return density;}
  double GetMolMass(){ return molarmass;}
  double GetAtomcm3(){ return atom_cm3;}//not quite atom_cm3 see the scale factor CHANGE????!!!
  double GetZ(){return z;}
  double GetMax(){return photoenergy.back();}
  double GetMin(){return photoenergy.front();}
  double GetBetheBloch(double);
  
  void GetImgDielectric();
  void GetRealDielectric();
  void GetRelCrossSection(double);

  bool GetReFlag(){ return set_real;}
  bool GetImgFlag(){ return set_img;}
  bool GetTableFlag(){ return set_table;}

  void SetPhotoEnergyVec(std::vector<double>);
  void SetPhotoValueVec(std::vector<double>);
  
  void SetReFlag(bool flag){ set_real = flag; }
  void SetImgFlag(bool flag){ set_img = flag; }
  void SetTableFlag(bool flag){ set_table = flag; }

  void WriteToFile(string);
  
  std::vector<double> GetEnergyVec(){ return photoenergy;}
  std::vector<double> GetPhotoCrossVec(){ return photovalue;}

  //  gsl_spline * GetPhotoCrossSpline(){return photo_cross_table;}
  //  Dielectric MixGas(Dielectric &, Dielectric &, double, double,int,double,double);
    Dielectric MixGas(Dielectric &, Dielectric &, double, double);

  double im_epsilon(double);
  double photo_cross_interp(double);
  double integprand(double);
  double real_interp(double);
  
  double InterpolateReal(double);
  double InterpolateImg(double);
  double InterpolateRelCrossSect(double);
  double integrand(double);
  double GetCrossSection(double,double,double);  
  double GetMoment(int,double);
  double GetRuthMoment(int);
  double GetRutherford(double);
  double GetConvolution(double,int);
  //  double dipole_oscill(double, void *);

  TGraph * DrawPhotoCross(int,double, double);
  TGraph * DrawImaginary();
  TGraph * DrawReal();
  TGraph * DrawCrossSection(bool);
  TGraph * DrawRutherford(int,double,double,bool);
  TGraph * DrawBichselSeg(double,double,double,double,double,double,int);
  TGraph * DrawConvolution(int,double,double,int,double);

  void ConvSelf(vector<double> &, vector<double> &);
  void ConvSelf(vector<double> &, vector<double> &, double);
  void ConvVec(vector<double> &, vector<double> &,vector<double> &, vector<double> &);
  double calcInt(vector<double> &, vector<double> &);
  TGraph * GraphFunc(vector<double> &, vector<double> &);
  
};
#endif  
