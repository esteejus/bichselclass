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
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *photo_cross_table;

  double energy_p = 0;
  double emax, emin;
  double avogadro    =  6.022e23;
  double hbar_c      = 1.97327e-5; // [eV *cm]
  double units_cross = 1e-18; //units the cross section is given in in cm^2
  double density = -1 , molarmass = -1 , atom_cm3 = -1;
  
  std::vector<double> photoenergy; //array to store values of cumulative dist
  std::vector<double> photovalue; //array to store values of cumulative dist

  std::vector<double> img_energy, img_value;
  std::vector<double> real_energy, real_value;

  bool set_table = false;
  bool set_img   = false;//if img is calculated or set through table this is TRUE
  bool set_real   = false;//if real is calculated or set through table this is TRUE
  
 public:
 Dielectric( double d_density, double d_molarmass) : density(d_density), molarmass(d_molarmass), atom_cm3( (d_density/d_molarmass) * (avogadro*units_cross) ){}
  
  void SetPhotoCross ( const std::string &);
  void SetDensity(double den){density = den;}
  void SetMolMass(double mm){molarmass = mm;}

  void GetImgDielectric(double,double,double);
  void GetRealDielectric(double,double,double);
  
  double im_epsilon(double);
  double photo_cross_interp(double);
  double integrand(double);
  //  double dipole_oscill(double, void *);

  TGraph * DrawPhotoCross(int,double, double);
  TGraph * DrawImaginary();
  TGraph * DrawReal();



};
#endif  
