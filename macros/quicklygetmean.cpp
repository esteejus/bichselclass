#include "bichsel.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;

vector<double> bgamma_vec,bgamma_value;

void SetBgammaTable (const std::string& filename) {
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      bgamma_vec.push_back(d_energy);
      bgamma_value.push_back(d_value);
    }
  }
  else std::cout << "Unable to open M0 data file" << std::endl;

  return;
}

int main(){

  SetBgammaTable("bgamma_pid.dat");

  gsl_interp_accel *acc1 = gsl_interp_accel_alloc();
  int size = bgamma_vec.size();
  gsl_spline *m0_table = gsl_spline_alloc(gsl_interp_akima,size);
  gsl_spline_init(m0_table,bgamma_vec.data(),bgamma_value.data(),size);
    
  double mass,mom,bgamma,z;
  cout<<"Enter mass "<<endl;
  cin>>mass;
  cout<<"Enter mom/z"<<endl;
  cin>>mom;
  cout<<"Enter z "<<endl;
  cin>>z;
  bgamma = (mom*z)/mass;
  cout<<bgamma<<endl;
  double value = gsl_spline_eval(m0_table,bgamma,acc1);
  cout<<"dE/dx is "<<endl;
  cout<<value<<endl;
  
  return 0;
}
