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


void SetTable(tk::spline &f_cross,const string& filename) {
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
      cout<<energy<<" "<<value<<endl;
    }
  }
  else cout << "File could not be opened" << endl;
  f_cross.set_points(x,y);

  return;
}


int main(){
  
  tk::spline ch4_real;
  tk::spline ar_real;
  tk::spline ar_blumrol_real;
  tk::spline p10_real;  
  
  SetTable(ch4_real,"./dielectricData/CH4_real.dat");
  //  SetTable(ar_real,"./dielectricData/Ar_real_qag.dat");
  //  SetTable(ar_real,"./dielectricData/Blum_rolandi_RealEpsilon.dat"); 
  SetTable(ar_real,"./dielectricData/Ar_real_west.dat"); 
  SetTable(p10_real,"./dielectricData/P10_real.dat");

  int num_points = 1e5;

  TGraph p10 = TGraph(num_points);
  TGraph ch4 = TGraph(num_points);
  TGraph ar = TGraph(num_points);
  TGraph ar_blumrol =fillgraph("./dielectricData/Blum_rolandi_RealEpsilon.dat");
  TGraph sum = TGraph(num_points);
  
  double max = 100.;
  double min = 8;
  double stepsize = (max-min)/num_points;
  
  double ar_f = .7;
  double ch4_f = .45;
  
  for(int i = 0;i < num_points; ++i)
    {
      double x = min + i*stepsize;
      p10.SetPoint(i,x,p10_real(x));
      ar.SetPoint(i,x,ar_real(x));
      ch4.SetPoint(i,x,ch4_real(x));
      sum.SetPoint(i,x,ar_real(x)*ar_f+ch4_real(x)*ch4_f);
    }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->cd();
  p10.GetXaxis()->SetRangeUser(8,100);
  p10.Draw();
  c1->SaveAs("./png/p10.png");
  
  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->cd();
  c2->SetLogx();
  ar_blumrol.GetXaxis()->SetRangeUser(8,100);
  ar_blumrol.Draw("APO");
  ar.SetMarkerColor(2);
  ar.Draw("same PO");
  c2->SaveAs("./png/ar.png");

  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->cd();
  ch4.GetXaxis()->SetRangeUser(8,100);
  ch4.Draw();
  c3->SaveAs("./png/ch4.png");

  TCanvas *c4 = new TCanvas("c4","c4",1);
  c4->cd();
  p10.GetXaxis()->SetRangeUser(8,100);
  p10.Draw();
  sum.SetMarkerColor(2);
  sum.Draw("same PO");
  c4->SaveAs("./png/sum_compare.png");

  return 0;

}
