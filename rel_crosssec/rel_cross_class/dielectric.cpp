#include "dielectric.h"

void Dielectric::SetPhotoCross (const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      photoenergy.push_back(d_energy);
      photovalue.push_back(d_value);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  double *x_a = photoenergy.data();
  double *y_a = photovalue.data();

  int size_a = photoenergy.size();
  emin = photoenergy.front();
  emax = photoenergy.back();

  photo_cross_table = gsl_spline_alloc(gsl_interp_akima,size_a);
  gsl_spline_init(photo_cross_table,x_a,y_a,size_a);

  //delete x_a;
  //  delete y_a;
  
  set_table = true;

  return;
}


void Dielectric::GetImgDielectric( double npoints, double x1=0, double x2=1.e3 ){
  double energy_step = (x2-x1)/npoints;

  for(int iEnergy = 1; iEnergy <= npoints; ++iEnergy){
    double energy = energy_step*iEnergy + x1;
    double value = im_epsilon(energy);
    img_energy.push_back(energy);
    img_value.push_back(value);
  }
  set_img = true;

  return;
}

double Dielectric::im_epsilon(double x) {
  if( (atom_cm3 + 1)<.0001 ){
    cout << "ERROR::Atomic density was not set properly set atom/cm3";
    return 0;
  }
  double coeff = (hbar_c * atom_cm3);

  return  (coeff*photo_cross_interp(x)/x) ;
}

double Dielectric::photo_cross_interp(double x) {
  double f=0.;//function value
  if ( set_table == false )
    {
      cout<<"Cross section table not set."<<endl;
      cout<<"Use SetPhotoCross() member function to set"<<endl;
      return -1;
    }

 if( x>=emin && x<=emax ) f =  gsl_spline_eval(photo_cross_table,x,acc);
 else f = 0;
 
  return f;
}

TGraph * Dielectric::DrawPhotoCross(int npoints, double x1=0, double x2=1.e3){
  double energy_step = (x2-x1)/npoints;
  TGraph * photocross = new TGraph(npoints);

  for(int iEnergy = 0 ; iEnergy < npoints; ++iEnergy)
    {
      double energy = energy_step*iEnergy + x1;
      double value = photo_cross_interp(energy);
      photocross -> SetPoint(iEnergy,energy,value);
    }

  return photocross;
  
}

TGraph * Dielectric::DrawImaginary(){
  int npoints = img_energy.size();
  TGraph * imaginary = new TGraph(npoints);

  //A couple of checks 
  if(img_energy.size() != img_value.size()){
    cout << "ERROR energy table and value table of the Imaginary componet arrays are not equal. Please check the file you are uploading as your table" <<endl;
    set_img = false;
}
  if( set_img == false ){
    cout<<"Imaginary values were not calculated or the tables were not set. Please use GetImgDielectric(double npoints) to calculate or set using ..."<<endl;
      return imaginary;
  }

  int img_size = img_energy.size();
for(int iImg = 0 ; iImg < img_size; ++iImg)
    imaginary -> SetPoint(iImg,img_energy.at(iImg),img_value.at(iImg));

  return imaginary;
  
}
