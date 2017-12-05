#include <iostream>

using namespace std;

int main(){

  double amu=931.5;     //MeV/c^2 from Nuclear Wallet cards
  //masses from nuclear wallet cards
  double p_mass = amu+7.289;
  double d_mass = amu*2+13.136;
  double t_mass = amu*3+14.95;
  double he3_mass = amu*3+14.931;
  double he4_mass = amu*4+2.425;

  
  cout<<"Proton mass "<<p_mass<<endl;
  cout<<"d mass "<<d_mass<<endl;
    cout<<"t mass "<<t_mass<<endl;
      cout<<"he3 mass "<<he3_mass<<endl;
      cout<<"he4 mass "<<he4_mass<<endl;
  return 0;
}
