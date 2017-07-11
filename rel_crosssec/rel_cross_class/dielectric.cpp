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

  //  double *x_a = photoenergy.data();
  //  double *y_a = photovalue.data();

  //  int size_a = photoenergy.size();
  //  emin = photoenergy.front();
  //  emax = photoenergy.back();

  //  photo_cross_table = gsl_spline_alloc(gsl_interp_akima,size_a);
  //  gsl_spline_init(photo_cross_table,photoenergy.data(),photovalue.data(),size_a);

  //delete x_a;
  //  delete y_a;
  
  set_table = true;

  return;
}

void Dielectric::SetPhotoEnergyVec(std::vector<double> vec){
	photoenergy.clear();
	photoenergy = vec;
	cout<<"Photo cross section energy values overwritten"<<endl;

    return;
  }

void Dielectric::SetPhotoValueVec(std::vector<double> vec){
	photovalue.clear();
	photovalue = vec;
	cout<<"Photo cross section values overwritten"<<endl;

    return;
  }

double Dielectric::GetMoment(int mom){
  int org_mom = moment;
  moment = mom;
  double result, error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::InterpolateRelCrossSect, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  gsl_integration_qags (F, 1, 1e3, 0, 1e-3, 1000,w, &result, &error);    
  gsl_integration_workspace_free (w);
  moment = org_mom;
  
  return ( atom_cm3 * z * result );
}

double Dielectric::GetRuthMoment(int mom){
  int org_mom = moment;
  moment = mom;
  double result, error;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::GetRutherford, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  gsl_integration_qags (F, 9, 1e4, 0, 1e-3, 1000,w, &result, &error);    
  gsl_integration_workspace_free (w);
  moment = org_mom;
  
  return (atom_cm3 * z * result);
}

double Dielectric::GetCrossSection(double x, double b, double intgrl){
  //b stands for beta as in beta*gamma
  double f = 0.;
  if(!(set_img == true && set_real == true))
    {
      cout<<"Real or img not calculated or set"<<endl;
      return -1;
    }
  double re_v        = InterpolateReal(x);
  double im_v        = InterpolateImg(x);
  double fine_struct  =  1./137;
  double coeff        = fine_struct/(pow(b,2) * 3.1415);
  double theta        = (im_v*pow(b,2))/(1-re_v*pow(b,2));// theta for phase
  double y_t          = (im_v*pow(b,2));// theta for phase
  double x_t          = (1-re_v*pow(b,2));// theta for phase
  //theta defined to be theta = atan2(y_t,x_t)
  double mag_eps      = pow(re_v,2)+pow(im_v,2);//actually magnitude squared
  double n_e          = atom_cm3 * z;//electron density
  //=============
  //first term
  //=============
  double log_term = pow( pow( 1-(pow(b,2)*re_v),2) + pow(b,4)*pow(im_v,2) ,-.5);
  f = (photo_cross_interp(x) * log(log_term))/(x*z);
  //======
  //Second term
  //======
  f += (photo_cross_interp(x)/(x*z))*log((2*m_elec*pow(b,2))/x);
  //=======
  //Third term
  //=======
  f += intgrl/(pow(x,2)*z);
  //Fourth term
  //======
  //f += (pow(b,2)-(re_v/mag_eps))*atan2(y_t,x_t)/(atom_cm3*hbar_c);
  f += (pow(b,2)-(re_v/mag_eps))*atan(theta)/(atom_cm3 * z * hbar_c);
  //  cout<<"tan "<<atan(theta)<<endl;
  //  cout<<"real "<<re_v<<endl;
  //  cout<<"first "<<(pow(b,2)-(re_v/mag_eps))<<endl;
    //cout<<"atom_cm3"<<atom_cm3<<endl;
    //  cout<<"f is "<<f<<endl;
  //=====
  //Multiply by overall coeff alpha/beta^2/pi
  //=====
  f *= coeff;
  

  return f;
  
}

void Dielectric::GetRelCrossSection(double b_gamma,int npoints, double x1, double x2){
  relcross_energy.clear();
  relcross_value.clear();

  double beta = b_gamma/sqrt(1+pow(b_gamma,2));
  double energy_step = (x2-x1)/npoints;
  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::photo_cross_interp, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  //  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  for(int iEnergy = 1; iEnergy <= npoints; ++iEnergy){
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double energy = energy_step*iEnergy + x1;
    double result, error;
    gsl_integration_qags (F, 0, energy, 0, 1e-3, 1000,w, &result, &error);    
    //Find the cross section
    double section = GetCrossSection(energy,beta,result);

    relcross_energy.push_back(energy);
    relcross_value.push_back(section);
    gsl_integration_workspace_free (w);
  }


  return;
}

void Dielectric::GetRealDielectric(double npoints, double x1=0, double x2=1.e3 ){
  real_energy.clear();
  real_value.clear();
  double energy_step = (x2-x1)/npoints;
  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::integrand, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
  for(int iEnergy = 1; iEnergy <= npoints; ++iEnergy){
    
    energy_p = energy_step*iEnergy + x1;
    double result, error;
    gsl_integration_qag (F, 0, 1e5, 1e-10, 1e-3,1000,6,w, &result, &error); 
    //      gsl_integration_qawc (&F, .1, 1e5, e, 0 , 1e-3, 1000,w, &result, &error); 
    real_energy.push_back(energy_p);
    real_value.push_back(result + 1);
  }
  gsl_integration_workspace_free (w);
  set_real = true;

  return;
}

void Dielectric::GetImgDielectric( double npoints, double x1=0, double x2=1.e3 ){
  img_energy.clear();
  img_value.clear();
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

double Dielectric::InterpolateReal(double energy){
  if ( set_real == false )
    {
      cout<<"Real dielectric not calculated or set."<<endl;
      return -1;
    }

  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  int size = real_energy.size();
  gsl_spline *real_table = gsl_spline_alloc(gsl_interp_akima,size);
  gsl_spline_init(real_table,real_energy.data(),real_value.data(),size);
  double value = gsl_spline_eval(real_table,energy,acc1);

  gsl_interp_accel_free(acc1);
  gsl_spline_free(real_table);
  
  return value;
}

double Dielectric::InterpolateImg(double energy){
  if ( set_img == false )
    {
      cout<<"ERROR: Imaginary dielectric not calculated or set."<<endl;
      return -1;
    }

  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  int size = img_energy.size();
  gsl_spline *img_table = gsl_spline_alloc(gsl_interp_akima,size);
  gsl_spline_init(img_table,img_energy.data(),img_value.data(),size);
  double value = gsl_spline_eval(img_table,energy,acc1);

  gsl_interp_accel_free(acc1);
  gsl_spline_free(img_table);
  
  return value;
}

//THIS IS NOTTTTTT PHOTO CROSS SECTION
double Dielectric::InterpolateRelCrossSect(double energy){
  //  if ( set_img == false )
  //    {
  //      cout<<"ERROR: Imaginary dielectric not calculated or set."<<endl;
  //      return -1;
  //    }

  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  int size = relcross_energy.size();
  gsl_spline *cross_table = gsl_spline_alloc(gsl_interp_akima,size);
  gsl_spline_init(cross_table,relcross_energy.data(),relcross_value.data(),size);
  double value = gsl_spline_eval(cross_table,energy,acc1);

  value *= pow(energy,moment);
  
  gsl_interp_accel_free(acc1);
  gsl_spline_free(cross_table);
  
  return value;
}

Dielectric Dielectric::MixGas(Dielectric &d1, Dielectric &d2, double d1_f, double d2_f, int npoints, double x1, double x2){

  double d3_density   = d1.GetDensity()*d1_f + d2.GetDensity()*d2_f;
  double d3_molarmass = d1.GetMolMass()*d1_f + d2.GetMolMass()*d2_f;
  double d3_z         = d1.GetZ()*d1_f + d2.GetZ()*d2_f;
  Dielectric d3(d3_density,d3_molarmass,d3_z);
  std::vector<double> d3_photoenergy, d3_photovalue;

  if(d1.GetTableFlag() == true && d2.GetTableFlag() == true)
    {
      double energy_step = (x2-x1)/npoints;
      for(int iEnergy = 0; iEnergy < npoints; ++iEnergy)
	{
	  double energy = energy_step*iEnergy + x1;
	  double d1_value = d1.photo_cross_interp(energy);
	  double d2_value = d2.photo_cross_interp(energy);

	  double real_mix = d1_value*d1_f + d2_value*d2_f;
	  cout<<"Energy in mix is "<<energy<<endl;
	  cout<<"D1, d2"<<d1_value<<" "<<d2_value<<endl;
	  d3_photoenergy.push_back(energy);
	  d3_photovalue.push_back(real_mix);
	}
      cout << "Photo absorbtion cross section of mixed gas set"<<endl;
    }
  else
    cout<<"WARNING::Could not set the mixed gas Photo absorbtion cross section. Check to make sure all gas componets tables  are set"<<endl;

  d3.SetPhotoEnergyVec(d3_photoenergy);  
  d3.SetPhotoValueVec(d3_photovalue);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  int size = d3_photoenergy.size();
  gsl_spline *photo_table = gsl_spline_alloc(gsl_interp_akima,size);
  gsl_spline_init(photo_table,d3_photoenergy.data(),d3_photovalue.data(),size);
  //  cout<<"Testing eval"<<gsl_spline_eval(photo_table,15,acc)<<endl;

      gsl_interp_accel_free(acc);
      gsl_spline_free(photo_table);

      d3.SetTableFlag(true);
      
  return d3;
}

double Dielectric::integrand(double x) {
  //def of QAWC integration in gsl library is
  //I = \int_a^b dx f(x) / (x - c)
  //we want to integrate I = \int_(0,infinity) dE' f(E') / (E'^2 - E^2)
  //which can also be simplified as f(E')/(E'-E)*(E'+E)
  //since we are integrating from 0 to infinity in E'
  //we can put f(E')/(E'+E)=f"(E')
  //thus the new integral is I = dE' f"(E')/(E'-E)
  //which is the same as the gsl format
  //with the singularity at (E'-E)

  //  struct f_params * params = (struct f_params *)p;
  double e = energy_p;
  
  double f = (x*im_epsilon(x)-e*im_epsilon(e));
  f = f/(pow(x,2)-pow(e,2));
  f = f*(2/3.1415);
    
  return f;
}

double Dielectric::im_epsilon(double x) {
  if( (atom_cm3 + 1)<.0001 ){
    cout << "ERROR::Atomic density was not set properly set atom/cm3";
    return 0;
  }
  double coeff = (hbar_c * atom_cm3);

  return  (coeff*photo_cross_interp(x)/x) ;
}

double Dielectric::real_interp(double x){
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  int size_a = real_energy.size();
  gsl_spline *real_table = gsl_spline_alloc(gsl_interp_akima,size_a);
  gsl_spline_init(real_table,real_energy.data(),real_value.data(),size_a);
  double f = gsl_spline_eval(real_table,x,acc);
  gsl_interp_accel_free(acc);
  gsl_spline_free(real_table);

  return f;

}
double Dielectric::photo_cross_interp(double x) {
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  int size_a = photoenergy.size();
  gsl_spline * photo_cross_table = gsl_spline_alloc(gsl_interp_akima,size_a);
  gsl_spline_init(photo_cross_table,photoenergy.data(),photovalue.data(),size_a);

  double f=0.;//function value
  if ( set_table == false )
    {
      cout<<"Cross section table not set."<<endl;
      cout<<"Use SetPhotoCross() member function to set"<<endl;
      return -1;
    }

  double emin = photoenergy.front();
  double emax = photoenergy.back();

  if( x>=emin && x<=emax ){ f =  gsl_spline_eval(photo_cross_table,x,acc);
  }
  else f = 0;

 gsl_interp_accel_free(acc);
 gsl_spline_free(photo_cross_table);

 return (f*units_cross);
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


TGraph * Dielectric::DrawReal(){
  int npoints = img_energy.size();
  TGraph * real = new TGraph(npoints);

  //A couple of checks 
  if(real_energy.size() != real_value.size()){
    cout << "ERROR energy table and value table of the Real componet arrays are not equal. Please check the file you are uploading as your table" <<endl;
    set_real = false;
}
  if( set_real == false ){
    cout<<"Real values were not calculated or the tables were not set. Please use GetImgDielectric(double npoints) to calculate or set using ..."<<endl;
      return real;
  }

  int real_size = real_energy.size();
for(int iRe = 0 ; iRe < real_size; ++iRe)
  real -> SetPoint(iRe,real_energy.at(iRe),real_value.at(iRe));

  return real;
  
}

TGraph * Dielectric::DrawCrossSection(bool mbarn){
  int npoints = relcross_energy.size();
  TGraph * cross = new TGraph(npoints);

for(int iCross = 0 ; iCross < npoints; ++iCross)
  {
    double value = relcross_value.at(iCross);
    if(mbarn == true)
      value *= 1e18;//show in Mb units

    cross -> SetPoint(iCross,relcross_energy.at(iCross),value);
  }


  return cross;
  
}

double Dielectric::GetRutherford(double energy){
  double bgamma = 3.6;
  double coeff = 2.5496e-19;//ev*cm^ 2
  //  double coeff = 1.5354e5;//ev*cm^ 2
  double beta = bgamma/sqrt(1+pow(bgamma,2));
  double emax = 2*m_elec*pow(bgamma,2);
  double value = coeff*(1-(pow(beta,2)*energy/emax))/pow(beta*energy,2);
  value *= pow(energy,moment);

  return value;
}

TGraph * Dielectric::DrawRutherford( int npoints, double x1, double x2,bool mbarn){
  double energy_step = (x2-x1)/npoints;
  TGraph * ruth = new TGraph(npoints);
  for(int iEnergy = 0 ; iEnergy < npoints; ++iEnergy)
    {
      double energy = energy_step*iEnergy + x1;
      double value = GetRutherford(energy);
      if(mbarn == true)
	value *= 1e18;//show in Mb units

      ruth -> SetPoint(iEnergy,energy,value);
    }
  
  return ruth;
}

