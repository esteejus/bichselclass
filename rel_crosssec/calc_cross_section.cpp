#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <gsl/gsl_integration.h>
#include "spline.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;


struct f_params {tk::spline b; double real; double imag; double beta; double den; double integrand;};

void SetTable(tk::spline &f_cross,std::vector<double> &x, std::vector<double> &y,const string& filename) {

  ifstream file (filename.c_str());
  double value = 0; // componet of dielectric function
  double energy = 0;
  string line;

  if (file.is_open()){
    while(getline (file,line) ){
      istringstream in(line);
      in>>energy;
      in>>value;
      x.push_back(energy);
      y.push_back(value);
    }
  }
  else cout << "File could not be opened" << endl;
  f_cross.set_points(x,y);

  return;
}

double photo_cross(double x, void * p) {
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);
  //if you want to make a nice code put the vectors used to interpolate
  //in to the first if statement vector.front() and .back() for the
  //piecewise function
  double f=0.;//function value
  if( x>=10.6406 && x<=500 ) f = f_cross(x);
  else if( x > 500 && x < 3206 ) f = 1.346*pow(500/x,2.54);
  else if( x >= 3206) f = .1*pow(3206/x,2.75);
  else f = 0;

  return f;

}

double cross_section(double x, void *p){
  double f = 0.;
  struct f_params * params = (struct f_params *)p;
  tk::spline f_cross = (params->b);
  double re_v        = (params->real);
  double im_v        = (params->imag);
  double b           = (params->beta);
  double density     = (params->den);
  double intgrl      = (params->integrand);
  double fine_struct  =  1./137;
  double coeff        = fine_struct/(pow(b,2) * 3.1415);
  double z            = 19;           //atomic number
  double m_elec       = 5.11e5;       //eV/c^2
  double hbar_c       = 1.9732697e-5; //eV*cm
  double theta        = (im_v*pow(b,2))/(1-re_v*pow(b,2));// theta for phase

  //=============
  //first term
  //=============
  double log_term = pow( pow(1-(pow(b,2)*re_v),2)+pow(pow(b,2)*im_v,2) ,-.5);
  f = (photo_cross(x,p) * log(log_term))/(x*z);
  //======
  //Second term
  //======
    f += (pow(b,2)-(re_v/(pow(re_v,2)+pow(im_v,2))))*atan(theta)/(density*hbar_c);
   //======
  //Third term
  //======
    f += (photo_cross(x,p)/(x*z))*log((2*m_elec*pow(b,2))/x);
  //=======
  //Fourth term
  //=======
  f += intgrl/(pow(x,2)*z);

  //=====
  //Multiply by overall coeff alpha/beta^2/pi
  //=====
  f *= coeff;
  

  return f;
  
}


int main(){
  //===================
  //GAS parameters
  //===================
  double density = 1.662e-3;   // [g/cm^3]
  double molarmass = 39.948;   // [g/mol]
  double avogadro = 6.022e23; // [atoms/mol]*scalefactor
  double atom_cm3 = (density / molarmass) * avogadro; //[atoms/cm^3] or N in the equation

  //=====================
  //Set the dilectric data files here
  //====================
  std::vector<double> e_re,e_im,re,im; //e_im and e_re are the energy values
  tk::spline imaginary,real;
  SetTable(imaginary,e_im,im,"imaginary.dat"); 
  SetTable(real,e_re,re,"real.dat"); 

  //=================
  //Set the photoelectric cross section
  //=================
  std::vector<double> e_fcross,fcross; //
  tk::spline f_cross;
  SetTable(f_cross,e_fcross,fcross,"./argon10eV_500eV.dat");

  //===================
  //Section for defining beta gamma
  //===================
  double b_gamma = 3.6;
  double beta = sqrt( pow(b_gamma,2)/(1+pow(b_gamma,2)) );

  //===================
  //Begin to find the reletavisitc cross seciton
  //==================
  int num_steps = 1e4;
  double e_stepsize = 1.e3/num_steps;
  cout<<"step size is "<<e_stepsize<<endl;
  
  std::vector<double>e_relcross,rel_cross;
  TGraph *g_cross = new TGraph(num_steps);

  for(int i=0;i<=num_steps;++i){
        double energy = i * e_stepsize;
	if(energy<10.7)continue;
cout<<"on "<<i<< "energy "<<energy<<endl;
    if(energy <= e_im.back() && energy <= e_re.back()){
      double ep_re = real(energy)*1.e-4 + 1; //table is Re[e]-1
      double ep_im = imaginary(energy)*1.e-4;
      
      //integrand is 0 because we are going to solve that in the next lines of code 
      struct f_params alpha = {f_cross,ep_re,ep_im,atom_cm3,beta,0};

      double result, error;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
      gsl_function F;
      F.function = &photo_cross;
      F.params = &alpha;
      gsl_integration_qags (&F, 10.6, energy, 0, 1e-3, 1000,w, &result, &error);    
      printf ("result          = % .18f\n", result);
      //      printf ("exact result    = % .18f\n", expected);
      printf ("estimated error = % .18f\n", error);
      //      printf ("actual error    = % .18f\n", result - expected);
      printf ("intervals       = %zu\n", w->size);

      gsl_integration_workspace_free (w);
      cout<<"result "<<result<<endl;

      //Find the cross section
      alpha = {f_cross,ep_re,ep_im,beta,atom_cm3,result}; //fill the value of the integrand
      void * p = &alpha;
      double section = cross_section(energy,p);
      e_relcross.push_back(energy);
      rel_cross.push_back(section);
      g_cross->SetPoint(i,energy,section);
      cout<<"section "<<section<<endl;    
      
    }
    else cout << "The energy value is out of bounds of the dielectric tables" <<endl;


  }

  //===========================
  //Integrating the cross section function to get M0 moment
  //===========================
    tk::spline s_cross;
  s_cross.set_points(e_relcross,rel_cross);

  TCanvas *c1 = new TCanvas(1);
  c1->SetLogx();
  c1->SetLogy();
  g_cross->GetXaxis()->SetRangeUser(10,500);
  g_cross->Draw();

  c1->SaveAs("rel_cross_section.jpg");
  
return 0;

}
