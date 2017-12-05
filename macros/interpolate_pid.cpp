#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "spline.h"

using namespace std;

void SetTable(tk::spline &f_cross,std::vector<double> &x, std::vector<double> &y,const string& filename) {

  ifstream file (filename.c_str());
  double d_x=0,d_y=0;
  string line;

  if (file.is_open()){
    while(getline (file,line) ){
      istringstream in(line);
      in>>d_x;
      in>>d_y;
      x.push_back(d_x);
      y.push_back(d_y);
    }
  }
  else cout << "File could not be opened" << endl;
  f_cross.set_points(x,y);

  return;
}

int main(){

  vector<double> betag,dedx,sigma;
  tk::spline energyloss, error;

  SetTable(energyloss,betag,dedx,"dedx.dat");
  //  SetTable(bg,sigma,error,"sigma.dat");

  ofstream out;
  out.open("t.dat");

  double mass = 3*938.;
  double start_m = 226;
  double end_m = 3000;
  int large_step = 40;
  int small_step = 30;
  double mom = start_m;

  while(mom<=end_m){

    double bg =0;
    double mom_atpoint3 = .3*mass;
    double mom_step_small = (mom_atpoint3 - start_m)/small_step;
    double mom_step_large = (end_m-mom_atpoint3)/large_step;

    if(mom<=mom_atpoint3) mom += mom_step_small;
    else mom += mom_step_large;

    bg = mom/mass;
    out<<mom<<"\t"<<energyloss(bg)<<endl;
  }

  return 0;

}
