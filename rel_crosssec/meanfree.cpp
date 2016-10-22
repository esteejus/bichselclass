#include <iostream>
#include <cmath>

using namespace std;

vector< vector<double> > SimpsonNInt(vector< vector <double> > &hist,double step_size,double low_lim, double up_lim){

  int num_steps = ceil((up_lim-low_lim)/step_size);//round up value

  std::vector<std::vector<double> > integral;
  //Initialize an array integral[i][j] where i runs from 0,1
  //0 componet stores the  x value in the cumulative integral [keV/cm]
  //1 componet stores the integral value which is normalized
  integral.resize(2);

  if(hist[0].size()!=num_steps)std::cout << "ERROR histogram and step size are different" << std::endl;

  if(hist[0].size()%2==0){std::cout<<"even"<<std::endl;hist[0].pop_back(); hist[1].pop_back();}//simpson's rule requies odd number set

  double sum = 0.;//cumulative integral
  double x = hist[0][0];//cumulative x posiiton
  
  for(int i=0;i<num_steps;i+=2){
    double dx =0.;//cumulative x position
    double df =0.;//integral step
    if( (i+2) != num_steps ){
      df = (step_size/3)*(hist[1][i] + 4*hist[1][i+1] + hist[1][i+2]);
      dx = (hist[0][i+2]-hist[0][i]);
      
      sum+=df;
      x+=dx;
      
      integral[0].push_back(x);
      integral[1].push_back(sum);
    }
    
    else break;
  }

  //last point is not a real physical point
  integral[0].pop_back();
  integral[1].pop_back();

  //scale the cumulative integral to 1
  //  for(int i=0;i<integral[0].size();++i) integral[1].at(i)=integral[1].at(i)/integral[1].back();

  return integral;
}



double rutherford(double e, double vel){
  double A = 40, Z =1 ;
  double cross_ruth = (

  return cross_ruth;
}

int main(){









  return 0;
}
