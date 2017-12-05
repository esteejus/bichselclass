#include "dielectric.h"

bool wayToSort(double i, double j) {return i < j;};

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
  set_table = true;

  return;
}

bool Dielectric::SetRealTable (const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open())
    {
      while(getline (file,line) ){
	std::istringstream in(line);
	in>>d_energy;
	in>>d_value;
	real_energy.push_back(d_energy);
	real_value.push_back(d_value);
      }
    }
  else
    {
      std::cout << "Unable to open real table file" << std::endl;
      return false;      
    }
  set_real = true;
  cout<<"Read in real values"<<endl;

  return true;;
}

bool Dielectric::SetImgTable (const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open())
    {
      while(getline (file,line) ){
	std::istringstream in(line);
	in>>d_energy;
	in>>d_value;
	img_energy.push_back(d_energy);
	img_value.push_back(d_value);
      }
    }
  else
    {
      std::cout << "Unable to open imaginary table file" << std::endl;
      return false;      
    }
  set_img = true;
  cout<<"Read in imaginary values"<<endl;  
  return true;;
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

double Dielectric::GetMoment(int mom,double bgamma){
  if(set_rel == false)
    {
      cout<<"Relatavistic cross section not calculated"<<endl;
      return -9999;
    }

  int org_mom = moment;
  moment = mom;
  double result, error;
  double emax = 2*m_elec*pow(bgamma,2);
  double beta = bgamma/sqrt(1+pow(bgamma,2));
  //  cout<<"Emax is "<<emax<<" beta is "<<beta<<endl;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::InterpolateRelCrossSect, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  gsl_integration_qag (F, 0, 1e6, 0, 1e-3, 1000,1,w, &result, &error);    
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
  gsl_integration_qag (F, 9, 1e4, 0, 1e-3, 1000,1,w, &result, &error);    
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
  //  double y_t          = (im_v*pow(b,2));// theta for phase
  //  double x_t          = (1-re_v*pow(b,2));// theta for phase
  //theta defined to be theta = atan2(y_t,x_t)
  double mag_eps      = pow(re_v,2)+pow(im_v,2);//actually magnitude squared
  //  double n_e          = atom_cm3 * z;//electron density
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
  //=====
  //Multiply by overall coeff alpha/beta^2/pi
  //=====
  f *= coeff;
  

  return f;
  
}

void Dielectric::GetRelCrossSection(double b_gamma){
  relcross_energy.clear();
  relcross_value.clear();

  double beta = b_gamma/sqrt(1+pow(b_gamma,2));
  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::photo_cross_interp, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  int npoints = photoenergy.size();
  for(int iEnergy = 0 ; iEnergy < npoints; ++iEnergy){
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double energy = photoenergy.at(iEnergy);
    double result, error;
    gsl_integration_qags (F, 0, energy, 0, 1e-3, 1000,w, &result, &error);    
    //Find the cross section
    double section = GetCrossSection(energy,beta,result);
    relcross_energy.push_back(energy);
    relcross_value.push_back(section);

    gsl_integration_workspace_free (w);
  }
  set_rel = true;

  return;
}

void Dielectric::GetRealDielectric(){
  real_energy.clear();
  real_value.clear();

  //Wrapper for the member function
  gsl_function_pp Fp( std::bind(&Dielectric::integrand, &(*this),  std::placeholders::_1) );
  gsl_function *F = static_cast<gsl_function*>(&Fp); 
  //  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  int npoints = photoenergy.size();  
  for(int iEnergy = 0; iEnergy < npoints; ++iEnergy){
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    energy_p = photoenergy.at(iEnergy)+.1;//due to endpoint issues add small shift like .1
    double result, error;
    gsl_integration_qag (F, 0, 1e5, 1e-10, 1e-2,1000,6,w, &result, &error); 
    real_energy.push_back(energy_p);
    real_value.push_back(result + 1);
    gsl_integration_workspace_free (w);
  }
  //  gsl_integration_workspace_free (w);
  set_real = true;
  cout<<"Real dielectric calculated"<<endl;

  return;
}

void Dielectric::GetImgDielectric(){
  img_energy.clear();
  img_value.clear();
  int npoints = photoenergy.size();
  for(int iEnergy = 0; iEnergy < npoints; ++iEnergy){
    double energy = photoenergy.at(iEnergy);
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
  double value = 0;
  if(energy >= real_energy.front() && energy <= real_energy.back())
    value = gsl_spline_eval(real_table,energy,acc1);
  else
    value = 1;//careful. Need better criteria here. Short fix
  
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
  double value = 0;
  if(energy >= img_energy.front() && energy <= img_energy.back())
    value = gsl_spline_eval(img_table,energy,acc1);
  else
    value = 0;//careful. Need better criteria here. Short fix

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
  double value = 0;
  if(energy >= relcross_energy.front() && energy <= relcross_energy.back())
    value = gsl_spline_eval(cross_table,energy,acc1);
  else
    value = 0;
  
  value *= pow(energy,moment);
  
  gsl_interp_accel_free(acc1);
  gsl_spline_free(cross_table);
  
  return value;
}

//Dielectric Dielectric::MixGas(Dielectric &d1, Dielectric &d2, double d1_f, double d2_f, int npoints, double x1, double x2){
Dielectric Dielectric::MixGas(Dielectric &d1, Dielectric &d2, double d1_f, double d2_f){
  double d3_density   = d1.GetDensity()*d1_f + d2.GetDensity()*d2_f;
  double d3_molarmass = d1.GetMolMass()*d1_f + d2.GetMolMass()*d2_f;
  double d3_z         = d1.GetZ()*d1_f + d2.GetZ()*d2_f;
  stringstream d3_name;
  d3_name<<d1.GetName()<<"_"<<d1_f<<d2.GetName()<<"_"<<d2_f;
  string str = d3_name.str();
  Dielectric d3(d3_density,d3_molarmass,d3_z,str);

  cout<<d3_density<<" "<<d3_molarmass<<" "<<d3_z<<" "<<str<<endl;
  std::vector<double> d3_photovalue;

  std::vector<double> d1_energy = d1.GetEnergyVec();
  std::vector<double> d2_energy = d2.GetEnergyVec();
  std::vector<double> d3_photoenergy(d1_energy.size() + d2_energy.size());
  std::vector<double>::iterator it;
  it=std::set_union (d1_energy.begin(), d1_energy.begin()+d1_energy.size(), d2_energy.begin(), d2_energy.begin()+d2_energy.size(), d3_photoenergy.begin());
  d3_photoenergy.resize(it-d3_photoenergy.begin());

  if(d1.GetTableFlag() == true && d2.GetTableFlag() == true)
    {
      //    double energy_step = (x2-x1)/npoints;
      int npoints = d3_photoenergy.size();
      for(int iEnergy = 0; iEnergy < npoints; ++iEnergy)
	{
	  double energy = d3_photoenergy.at(iEnergy);
	  //Divide by scale factor because photo_cross_interp accepts units of Mb
	  double d1_value = d1.photo_cross_interp(energy)/units_cross;
	  double d2_value = d2.photo_cross_interp(energy)/units_cross;

	  double real_mix = d1_value*d1_f + d2_value*d2_f;
	  d3_photovalue.push_back(real_mix);
	}
      cout << "Photo absorbtion cross section of mixed gas set"<<endl;
    }
  else
    cout<<"WARNING::Could not set the mixed gas Photo absorbtion cross section. Check to make sure all gas componets tables  are set"<<endl;

  d3.SetPhotoEnergyVec(d3_photoenergy);  
  d3.SetPhotoValueVec(d3_photovalue);
  /*
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  int size = d3_photoenergy.size();
  gsl_spline *photo_table = gsl_spline_alloc(gsl_interp_akima,size);
  gsl_spline_init(photo_table,d3_photoenergy.data(),d3_photovalue.data(),size);

  gsl_interp_accel_free(acc);
  gsl_spline_free(photo_table);
  */
  
  d3.SetTableFlag(true);
  cout<<"Done mixing gas"<<endl;
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
  if( x>=emin && x<=emax )
    f =  gsl_spline_eval(photo_cross_table,x,acc);
  else
    f = 0;
  gsl_interp_accel_free(acc);
  gsl_spline_free(photo_cross_table);

  return (f*units_cross);
}


double Dielectric::GetBetheBloch(double bgamma){
  double beta = bgamma/sqrt(1+pow(bgamma,2));
  double ioniz = 180;
  double tmax = 2*m_elec*pow(bgamma,2);
  double bb =0;
  double coeff = .307075e6;//[eV*cm^2/g]
  coeff *= z/(molarmass*pow(beta,2));
  bb = (2*m_elec*pow(bgamma,2)*tmax)/pow(ioniz,2);
  bb = .5*log(bb)-pow(beta,2);
  bb *= coeff;
  
  return (bb*density/1.e3);
}

TGraph * Dielectric::DrawPhotoCross(int npoints, double x1=0, double x2=1.e3){
  double energy_step = (x2-x1)/npoints;
  TGraph * photocross = new TGraph(npoints);

  for(int iEnergy = 0 ; iEnergy < npoints; ++iEnergy)
    {
      double energy = energy_step*iEnergy + x1;
      double value = photo_cross_interp(energy);
      //      cout<<energy<<" "<<value<<endl;
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
  for(int iRe = 0 ; iRe < real_size; ++iRe){
    real -> SetPoint(iRe,real_energy.at(iRe),real_value.at(iRe));
  }
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
      //    value*=relcross_energy.at(iCross);//for displaying E*cross
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

void Dielectric::WriteToFile(string opt){
  ofstream output;

  if(opt == "real")
    {
      stringstream out;
      out<<"real_"<<name<<".dat";
      string str = out.str();
      output.open(str);
      int size = real_energy.size();
      for(int i = 0 ;i < size;++i)
	output<<real_energy.at(i)<<"\t"<<real_value.at(i)<<endl;
      cout<<str<<" Wrote to file"<<endl;
    }
  else if(opt == "img")
    {
      stringstream out;
      out<<"img_"<<name<<".dat";
      string str = out.str();
      output.open(str);
      int size = img_energy.size();
      for(int i = 0 ;i < size;++i)
	output<<img_energy.at(i)<<"\t"<<img_value.at(i)<<endl;
      cout<<str<<" Wrote to file"<<endl;
    }
  else if(opt == "photocross")
    {
      stringstream out;
      out<<"photocross_"<<name<<".dat";
      string str = out.str();
      output.open(str);
      int size = photoenergy.size();
      for(int i = 0 ;i < size;++i)
	output<<photoenergy.at(i)<<"\t"<<photovalue.at(i)<<endl;
      cout<<str<<" Wrote to file"<<endl;
    }
  else
    cout<<"Enter a valid option for WriteToFile"<<endl;
  

  return;
}

double Dielectric::GetConvolution(double delta,int n){

  if(n ==1)
    return InterpolateRelCrossSect(delta)/units_cross;
  else
    {
      auto fcn = [&](double e)->double{
	double sigma_e = InterpolateRelCrossSect(e)/units_cross;
	double sigma_d = GetConvolution(delta - e, n-1);
	//	if(n=3)	cout<<"value"<<sigma_e*sigma_d<<endl;
	return sigma_e*sigma_d;};

      double error,result;
      std::function<double(double)> F_int(fcn);
      gsl_function_pp F1(F_int);
      gsl_function *F= static_cast<gsl_function*>(&F1); 
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
      gsl_integration_qag (F, 0, delta, 0, 1e-2, 1000,1,w, &result, &error);    
      gsl_integration_workspace_free (w);
      
      return result;
    }

}

TGraph * Dielectric::DrawConvolution(int npoints, double x1, double x2, int n,double scale){

  TGraph *conv = new TGraph(npoints);
  double energy_step = (x2-x1)/npoints;

  for(int iEnergy = 0; iEnergy < npoints; iEnergy++)
    {
      double delta = iEnergy*energy_step + x1 ;
      double value = GetConvolution(delta,n)*scale;
      cout<<"n, delta "<<n<<" "<<delta<<" "<<value<<endl;
      conv -> SetPoint(iEnergy,delta,value);
    }
  
  return conv;
}

TGraph * Dielectric::GraphFunc(vector<double> &a, vector<double> &b){
  int size_x = a.size();
  int npoints = 2*a.size();
  TGraph * f = new TGraph(npoints);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
  
  double x1 = a.front(), x2 = a.back();
  if(x1 < 1e-3) x1 = .1;//using log cant have nan errors
  double log_energy_step = (log(x2)-log(x1))/npoints;
  for(int i = 0 ;i < npoints; i++)
    {
      double delta = i*log_energy_step + log(x1) ;
      delta = exp(delta);
      double value = 0;
      if(delta < a.back() && delta >a.front())
	value = gsl_spline_eval(dist_table,delta,acc);

      //      if(value<0)cout<<"Value interpolat is neg"<<delta<<" "<<value<<endl;
      f->SetPoint(i,delta,value);
    }
  
  return f;
}

void Dielectric::GetCDF(vector<double> &a, vector<double> &b,vector<double>&cdf_x, vector<double>&cdf_y){
  cdf_x.clear();
  cdf_y.clear();  

  int size_x = a.size();
  int npoints = 2*a.size();

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
      
  auto integral = [&](double x)->double{
    double f = 0;
    if(x >= a.front() && x <= a.back())
      f = gsl_spline_eval(dist_table,x,acc);
    else
      f = 0;

    return f;};

  double x1 = a.front(), x2 = a.back();
  if(x1 < 1e-3) x1 = .1;//using log cant have nan errors
  double log_energy_step = (log(x2)-log(x1))/npoints;
  double prev_result = 0;
  double result = 0,error = 0;
  for(int i = 0 ;i < npoints; i++)
    {
      double delta = i*log_energy_step + log(x1) ;
      delta = exp(delta);

      std::function<double(double)> F2(integral);
      gsl_function_pp F2_2(F2);
      gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
      gsl_integration_qag (F_2, x1, delta, 0, 1e-3, 10000,6,w, &result, &error);    
      gsl_integration_workspace_free (w);
      if(result > prev_result)
	{
	  cdf_x.push_back(delta);
	  cdf_y.push_back(result);
	  prev_result = result;
	}

    }
      gsl_interp_accel_free(acc);
      gsl_spline_free(dist_table);

      return;
}

/*
void Dielectric::GraphCDF(vector<double> &a, vector<double> &b,vector<double>&cdf_x, vector<double>cdf_y){
  int size_x = a.size();
  int npoints = 2*a.size();

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
      
  auto integral = [&](double x)->double{
    double f = 0;
    if(x >= a.front() && x <= a.back())
      f = gsl_spline_eval(dist_table,x,acc);
    else
      f = 0;

    return f;};

  TGraph * f = new TGraph(npoints);
  
  double x1 = a.front(), x2 = a.back();
  if(x1 < 1e-3) x1 = .1;//using log cant have nan errors
  double log_energy_step = (log(x2)-log(x1))/npoints;
  double prev_result = 0;
  double result = 0,error = 0;
  for(int i = 0 ;i < npoints; i++)
    {
      double delta = i*log_energy_step + log(x1) ;
      delta = exp(delta);

      std::function<double(double)> F2(integral);
      gsl_function_pp F2_2(F2);
      gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
      gsl_integration_qag (F_2, x1, delta, 0, 1e-2, 10000,6,w, &result, &error);    
      gsl_integration_workspace_free (w);
      if(result > prev_result)
	{
	  f->SetPoint(i,delta,result);
	  cout<<result<<" "<<endl;
	  prev_result = result;
	}
      else
	{
	  prev_result = result;
	  continue;
	}
    }
      gsl_interp_accel_free(acc);
      gsl_spline_free(dist_table);

  return f;
}
*/
double Dielectric::calcInt(vector<double> &a, vector<double> &b){
  double x1 = a.front();
  double x2 = a.back();  
  int size_x = a.size();
  double result = 0, error = 0;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
      
  auto integral = [&](double x)->double{
    double f = 0;
    if(x >= a.front() && x <= a.back())
      f = gsl_spline_eval(dist_table,x,acc);
    else
      f = 0;

    return f;};

  std::function<double(double)> F2(integral);
  gsl_function_pp F2_2(F2);
  gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  gsl_integration_qag (F_2, x1, x2, 0, 1e-2, 10000,3,w, &result, &error);    

  gsl_integration_workspace_free (w);
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);

  return result;
}

void Dielectric::ConvSelf(vector<double> &h_delta, vector<double> &h_value,double c1){
  int size_x = h_delta.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  gsl_spline_init(dist_table,h_delta.data(),h_value.data(),size_x);

  for(int iEnergy = 0 ;iEnergy < size_x; iEnergy++){
    double result = 0, error = 0;
    double delta = h_delta.at(iEnergy) ;

    //calculate Int{h(delta-x)h(x)}
    auto integral = [&](double x)->double{
      double f = 0, f_d = 0;
      if(x >= h_delta.front() && x <= h_delta.back())
	f = gsl_spline_eval(dist_table,x,acc);
      else
	f = 0;
      if((delta-x) >= h_delta.front() && (delta-x) <= h_delta.back())
	f_d = gsl_spline_eval(dist_table,delta-x,acc);
      else
	f_d = 0;
      //      if(f_d*f<0)cout<<"f "<<"f_d "<<f<<" "<<f_d<<" "<<x<<" "<<delta-x<<endl;
      return f*f_d;};

    std::function<double(double)> F2(integral);
    gsl_function_pp F2_2(F2);
    gsl_function *F_2= static_cast<gsl_function*>(&F2_2); 
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    gsl_integration_qag (F_2, 0, delta, 1e-10, 1e-2, 1000,3,w, &result, &error);    
    double h_d = gsl_spline_eval(dist_table,delta,acc);
    result += 2*c1*h_d;
    h_value.at(iEnergy)=result;
    gsl_integration_workspace_free (w);
  }
    gsl_interp_accel_free(acc);
    gsl_spline_free(dist_table);
      
}


TGraph * Dielectric::DrawBichselSeg(double mass, double mom, double seg, double npoints, double x1, double x2){

  cout<<"=============================="<<endl;
  cout<<"Begining to calculate Bichsel distribution"<<endl;
  cout<<"For segment of length "<<seg<<" [cm]"<<endl;
  cout<<"=============================="<<endl;

  double bgamma = mom/mass;
  //mass in MeV/c^2
  //segment in cm
  relcross_energy.clear();
  relcross_value.clear();
  GetRelCrossSection(bgamma);

  double m_0 = GetMoment(0,bgamma);
  double n_0 = .001;                              //want to get n_0 to be approx this. 
  int pow_2 = ceil( log(seg * m_0 / n_0)/log(2) );//rounding up makes sure n_0 is at least above value
  double dx = seg/pow(2,pow_2);                   //step size in segment
  n_0 =  dx * m_0;                                //update n_0 actual n_0 used is dx*m_0
  
  if(x1 < 1e-4)x1 = .1;                           //using log steps x1=0 not valid input
  double log_energy_step = (log(x2)-log(x1))/npoints;

  //starting initial distribution set up
  std::vector<double> dist_vec, dist_value;
  std::vector<double> int_vec, int_value;         //initial dist

  double c1 = (1-n_0);
  double x_t = dx;                                //current segment length in convolution
  for(int iEnergy = 0; iEnergy < npoints; iEnergy++)
    {
      double delta = iEnergy*log_energy_step + log(x1) ;
      delta = exp(delta);
      double result = InterpolateRelCrossSect(delta)*dx;
      result *= atom_cm3 * z;//equation 7.1 in Bichsel NIM A references sigma which is cross*z*atom_cm3
      dist_vec.push_back(delta);
      dist_value.push_back(result);
    }
  cout<<"End of initialization"<<endl;

  while( 2*x_t <= seg )//need some kind of less than delta here
    {
	  ConvSelf(dist_vec,dist_value,c1);
	  x_t *= 2;
	  c1 = pow(c1,2);
	  cout<<"x is at "<<x_t<<endl;
    }

  //  double scale = calcInt(dist_vec,dist_value);
  //  cout<<"Integral is "<<scale<<endl;
  /*
  for(int i = 0; i < dist_vec.size(); i++)
    {
      dist_vec.at(i) /= 1000;
      //   dist_value.at(i) /= dist_vec.at(i);
    }
  double scale = calcInt(dist_vec,dist_value);
  for(int i = 0; i < dist_vec.size(); i++)
    {
      dist_value.at(i) /= scale;
      //    dist_value.at(i) *= 2.19269;
    }
  */
  TGraph *dist = GraphFunc(dist_vec,dist_value);
  //  TGraph *dist = GraphCDF(dist_vec,dist_value);
  return dist;
}

TH1D * Dielectric::GetMCdist(TGraph *f_dist,TString name,double seg_length, int num_seg, double t_frac,int num_mc_events)
{
  //frac is the truncaiton fraction i.e. .7 would mean lowest 70% kept and highest 30% is not
  cout<<"Starting MC section"<<endl;
  vector<double> dist_x, dist_y;
  int npoints = f_dist->GetN();
  for(int iGraph = 0 ;iGraph < npoints; iGraph++)
    {
      double x = 0, y = 0;
      f_dist->GetPoint(iGraph,x,y);
      dist_x.push_back(x);
      dist_y.push_back(y);
    }

  vector<double> cdf_x, cdf_y;
  GetCDF(dist_x,dist_y,cdf_x,cdf_y);
  //  cout<<"IS IT MONOTONIC? "<<std::is_sorted(cdf_y.begin(),cdf_y.end())<<endl;

TRandom3 ran = TRandom3(0);//0 seed will be based on clock of comp and give diff seed each session

 double dE = 10.; //[eV]
 double emax = cdf_x.back()- 60;//go 6 bins (10eV each) under min
 double emin = cdf_x.front()+ 60;//6 bins over max
 double r_max = cdf_y.back();//random number max ~=1
 double r_min = cdf_y.front();//randome number min ~= 0

 int nbins = ceil( (emax-emin)/dE );
 emax = nbins*dE + emin;

TH1D * c_dist = new TH1D(name,name,nbins,emin,emax);

int size_cdf = cdf_y.size();
gsl_interp_accel *acc = gsl_interp_accel_alloc();
gsl_spline *cdf_table = gsl_spline_alloc(gsl_interp_akima,size_cdf);
gsl_spline_init(cdf_table,cdf_y.data(),cdf_x.data(),size_cdf);


 for(int iMC = 0; iMC < num_mc_events; iMC++)
   {
     vector<double> track;
     for(int iSeg = 0 ;iSeg < num_seg; iSeg++)
       {
	 double rand = ran.Uniform(r_min,r_max);
	 double energy_loss = 0;
	 energy_loss = gsl_spline_eval(cdf_table,rand,acc);
	 track.push_back(energy_loss/seg_length);
       }
     
     sort(track.begin(),track.end(),wayToSort);   //sort accending
     int elem = floor( t_frac * track.size() );   //truncate; element to truncate from is elem
     track.erase(track.begin()+elem, track.end());//erase upper fraction
     double sum = std::accumulate(track.begin(), track.end(), 0.0);
     double mean = sum/track.size();
     c_dist -> Fill(mean);
   }

 double scale = c_dist -> Integral("width");
 c_dist -> Scale(1./scale);

 cout<<"End of MC section"<<endl;
  return c_dist;

}


TGraph * Dielectric::GetCScaling(TH1 *h1, TH1 *h2){

  vector<double> h1_x, h1_y, h2_x, h2_y;
  
  int binh1 = h1->GetXaxis()->GetNbins();
  int binh2 = h2->GetXaxis()->GetNbins();  

  cout<<"IN scaling section "<<binh1<<" "<<binh2<<endl;
  double prev_cdf = -1;
  for(int i = 2 ;i < binh1; i++)
    {
      double h1_x_value = h1 -> GetBinCenter(i);
      double h1_cdf = h1 -> Integral(1,i,"width");
      if(h1_cdf > prev_cdf)//this ensures monotomic 
	{
	  h1_x.push_back(h1_x_value);
	  h1_y.push_back(h1_cdf);
	  prev_cdf = h1_cdf;
	}
    }
  prev_cdf = -1;
  for(int i = 2 ;i < binh2; i++)
    {
      double h2_x_value = h2 -> GetBinCenter(i);
      double h2_cdf = h2 -> Integral(1,i,"width");
      if(h2_cdf > prev_cdf)//this ensures monotomic 
	{
	  h2_x.push_back(h2_x_value);
	  h2_y.push_back(h2_cdf);
	  prev_cdf = h2_cdf;
	}
    }

  cout<<"Finished filling vectors"<<endl;
int size_h2 = h2_y.size();
gsl_interp_accel *acc_h1 = gsl_interp_accel_alloc();
gsl_spline *h2_table = gsl_spline_alloc(gsl_interp_akima,size_h2);
gsl_spline_init(h2_table,h2_y.data(),h2_x.data(),size_h2);

int size_h1 = h1_y.size();
gsl_interp_accel *acc_h2 = gsl_interp_accel_alloc();
gsl_spline *h1_table = gsl_spline_alloc(gsl_interp_akima,size_h1);
gsl_spline_init(h1_table,h1_y.data(),h1_x.data(),size_h1);
 

 int npoints = 1e3;
 TGraph *graph = new TGraph(npoints);
 double step_size = 1./npoints;
 for(int i = 0 ;i < npoints; i++)
   {
     double x = step_size*i;
     double energy_h1 = gsl_spline_eval(h1_table,x,acc_h1);
     double energy_h2 = gsl_spline_eval(h2_table,x,acc_h2);
     graph -> SetPoint(i,energy_h1,energy_h2);
   }

 //  TGraph *graph = new TGraph(h1_x.size(),h1_x.data(),h1_y.data());

  return graph;
}
