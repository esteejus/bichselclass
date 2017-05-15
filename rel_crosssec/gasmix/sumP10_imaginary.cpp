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


void SetTable(vector <double> &x, vector <double> &y,const string& filename) {
  //  std::vector<double> x,y;
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
      cout<<energy<<" "<<value<<endl;
    }
  }
  else cout << "File could not be opened" << endl;

  return;
}


int main(){
  
  vector <double> ch4_img_e;
  vector <double> ch4_img_v;

  vector <double> ar_img_e;
  vector <double> ar_img_v;

  SetTable(ch4_img_e,ch4_img_v,"./dielectricData/CH4_Img.dat");
  SetTable(ar_img_e,ar_img_v,"./dielectricData/Ar_Img.dat");

  int num_points = ar_img_e.size();
  ofstream output;
  output.open("P10_Img.dat");

  cout<<"NUM OF POINTS"<<num_points<<endl;
  
  TGraph ch4 = TGraph(num_points);
  TGraph ar = TGraph(num_points);
  TGraph sum = TGraph(num_points);

  double ar_f = .9;
  double ch4_f = .1;
  
  if(ar_img_e.size()!=ch4_img_e.size())cout<<"ARRAYS NOT EQUAL IN SIZE"<<endl;
  
  for(int i=0;i<ar_img_e.size();++i)
    {
      if(!(ar_img_e.at(i)-ch4_img_e.at(i)<.0001))cout<<"Energies in array not equal"<<endl;

      ar.SetPoint(i,ar_img_e.at(i),ar_img_v.at(i));
      ch4.SetPoint(i,ch4_img_e.at(i),ch4_img_v.at(i));
      double p10_img = ar_f*ar_img_v.at(i) + ch4_f*ch4_img_v.at(i);

      sum.SetPoint(i,ar_img_e.at(i),p10_img);
      output<<ar_img_e.at(i)<<"\t"<<p10_img<<endl;
    }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->cd();
  c1->SetLogy();
  c1->SetLogx();
  sum.GetXaxis()->SetRangeUser(1,1e4);
  sum.Draw();
  c1->SaveAs("./png/sum_img.png");
  
  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->cd();
  c2->SetLogy();
  c2->SetLogx();
  ar.GetXaxis()->SetRangeUser(1,1e4);
  ar.Draw();
  c2->SaveAs("./png/ar_img.png");

  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->cd();
  c3->SetLogy();
  c3->SetLogx();
  ch4.GetXaxis()->SetRangeUser(1,1e4);
  ch4.Draw();
  c3->SaveAs("./png/ch4_img.png");

  return 0;

}
