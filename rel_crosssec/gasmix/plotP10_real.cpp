#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "spline.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

gsl_interp_accel *acc = gsl_interp_accel_alloc();

TGraph fillgraph(const string& filename) {
  std::vector<double> x,y;
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

  int npoints = x.size();
  TGraph graph = TGraph(npoints);
  for(int i = 0;i < npoints;++i)graph.SetPoint(i,x.at(i),y.at(i));

  return graph;
}


gsl_spline* SetTable(const string& filename) {
  gsl_spline *f_cross = NULL;
  std::vector<double> x,y;
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
      //      cout<<energy<<" "<<value<<endl;
    }
  }
  else cout << "File could not be opened" << endl;

  double *x_a = x.data();
  double *y_a = y.data();

  int size_a = x.size();
  f_cross = gsl_spline_alloc(gsl_interp_akima,size_a);
  gsl_spline_init(f_cross,x_a,y_a,size_a);

  return f_cross;
}


int main(){

ofstream output;
 output.open("P10_real.dat");
  
  gsl_spline *ch4_real =  SetTable("./dielectricData/CH4_real.dat");
  gsl_spline *ar_real =  SetTable("./dielectricData/Ar_real_mc.dat"); 

  int num_points = 5e3;

  TGraph ch4 = TGraph(num_points);
  TGraph ar = TGraph(num_points);
  TGraph sum = TGraph(num_points);
  
  double max = 9.9e2;
  double min = 1.2;
  double stepsize = (max-min)/num_points;
  
  double ar_f = .9;
  double ch4_f = .1;
  
  for(int i = 0;i < num_points; ++i)
    {
      double x = min + i*stepsize;

      double ar_v = gsl_spline_eval(ar_real,x,acc);
      double ch4_v = gsl_spline_eval(ch4_real,x,acc);

      if(ar_v>0)cout<<"HERE IS ERROR"<<x<<" "<<ar_v<<endl;      
      ar.SetPoint(i,x,ar_v);
      ch4.SetPoint(i,x,ch4_v);
      double p10_sum = ar_f*ar_v+ch4_f*ch4_v;
      sum.SetPoint(i,x,p10_sum);

      output<<x<<"\t"<<p10_sum<<endl;
    }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->cd();
  //  c2->SetLogy();
  //  c2->SetLogx();
  sum.GetXaxis()->SetRangeUser(1,1e2);
  sum.Draw();
  c1->SaveAs("./png/sum_real_plot.png");

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->cd();
  //  c2->SetLogy();
  c2->SetLogx();
  ar.GetYaxis()->SetRangeUser(-.0005,.001);
  ar.GetXaxis()->SetRangeUser(1,1e3);
  ar.Draw();
  c2->SaveAs("./png/ar_real_plot.png");

  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->cd();
  //  c3->SetLogy();
  //  c3->SetLogx();
  ch4.GetXaxis()->SetRangeUser(1,1e2);
  ch4.Draw();
  c3->SaveAs("./png/ch4_real_plot.png");

  return 0;

}
