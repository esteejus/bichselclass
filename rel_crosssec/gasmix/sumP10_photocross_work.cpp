#include <gsl/gsl_integration.h>
//#include <gsl/gsl_erron.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "spline.h"
#include <cmath>
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
using namespace std;

gsl_interp_accel *acc = gsl_interp_accel_alloc();

void SetTable(gsl_spline *f_cross, vector<double> &x , vector<double> &y ,const string& filename) {
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
  
  return;
}


int main(){
  ofstream output;
  output.open("P10_photocross.dat");
  std::vector<double> ch4_x;
  std::vector<double> ch4_y;
  std::vector<double> ar_x,ar_y; //vectors to store the table data

 gsl_spline *ch4_cross;
 gsl_spline *ar_cross;

 SetTable(ch4_cross,ch4_x,ch4_y,"../ch4_fisyak.dat");
 SetTable(ar_cross,ar_x,ar_y,"../argon_10_500_Sakamoto.dat"); 

 double * ch_x = &ch4_x[0];
 double * ch_y = &ch4_y[0];

 double * ar_a = &ar_x[0];
 double * ar_b = &ar_y[0];

 int size_ch = ch4_x.size();
 int size_ar = ar_x.size();

  ch4_cross = gsl_spline_alloc(gsl_interp_akima,size_ch);
  gsl_spline_init(ch4_cross,ch_x,ch_y,size_ch);

  ar_cross = gsl_spline_alloc(gsl_interp_akima,size_ar);
  gsl_spline_init(ar_cross,ar_a,ar_b,size_ar);

  int num_points = 5e4;

  TGraph ch4 = TGraph(num_points);
  TGraph ar = TGraph(num_points);
  TGraph sum = TGraph(num_points);
  
  //  TGraph *ch4_t = fillgraph("../ch4_fisyak.dat");
  //  TGraph *ch4_s = fillgraph("../ch4.dat");
  //  TGraph *ar_p = fillgraph("../argon_10_500_Sakamoto.dat"); 

  double max = 9.e3;
  double min = 11;
  double stepsize = (max-min)/num_points;
  
  double ar_f = .9;
  double ch4_f = .1;

  for(int i = 0;i < num_points; ++i)
    {
      double x = min + i*stepsize;

      double ar_v = gsl_spline_eval(ar_cross,x,acc);
      double ch4_v = gsl_spline_eval(ch4_cross,x,acc);

      if(ar_v < 0 || ch4_v < 0)cout<<"HERE IS ERROR"<<x<<" "<<ch4_v<<endl;      
      ar.SetPoint(i,x,ar_v);
      ch4.SetPoint(i,x,ch4_v);
      double p10_sum = ar_f*ar_v+ch4_f*ch4_v;
      sum.SetPoint(i,x,p10_sum);

      output<<x<<"\t"<<p10_sum<<endl;
    }

  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->cd();
  c1->SetLogy();
  c1->SetLogx();
  //  sum.GetXaxis()->SetRangeUser(1,1e2);
  sum.Draw();
  c1->SaveAs("./png/p10_cross_plot.png");

  TCanvas *c2 = new TCanvas("c2","c2",1);
  c2->cd();
  c2->SetLogy();
  c2->SetLogx();
  ar.GetYaxis()->SetRangeUser(0,1e2);
  //  ar.GetXaxis()->SetRangeUser(1,1e3);
  ar.Draw();
  c2->SaveAs("./png/ar_cross_plot.png");
  
  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->cd();
  c3->SetLogy();
  c3->SetLogx();
  //ch4.GetXaxis()->SetRangeUser(.9e2,1e2);
  ch4.Draw();
  c3->SaveAs("./png/ch4_cross_plot.png");
  /*
TCanvas *c4 = new TCanvas("c4","c4",1);
  c4->cd();
  c4->SetLogy();
  c4->SetLogx();
  ch4_t->GetXaxis()->SetRangeUser(275,300);
  ch4_t->Draw("");
  ch4_s->SetMarkerColor(2);
  ch4_s->Draw("SAME PO");
  c4 -> SaveAs("./png/pointsCH4.png");

TCanvas *c5 = new TCanvas("c5","c5",1);
  c5->cd();
  c5->SetLogy();
  c5->SetLogx();
  //  ch4_t->GetXaxis()->SetRangeUser(275,300);
  ar_p->Draw("");
  c5 -> SaveAs("./png/pointsAR.png");
  */
  
  return 0;

}



