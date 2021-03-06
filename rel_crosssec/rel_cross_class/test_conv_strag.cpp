#include "dielectric.h"

using namespace std;

vector<double>photoenergy,photovalue;
vector<double>cdf_x, cdf_y;

TGraph * GraphCDF2Dist(vector<double> &a, vector<double> &b){
  int size_x = a.size();
  int npoints = 100*a.size();
  TGraph * f = new TGraph(npoints);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_cspline,size_x);
  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
  
  double x1 = a.front(), x2 = a.back();
  double energy_step = (x2 - x1)/npoints;
  cout<<"energy step "<<energy_step<<endl;
  for(int i = 0 ;i < npoints; i++)
    {
      double delta_1 = i*energy_step + x1 ;
      double delta_2 = (i+1)*energy_step + x1 ;
      double y1 = 0, y2 = 0;
      if(delta_1 < a.back() && delta_1 >a.front())
	y1 = gsl_spline_eval(dist_table,delta_1,acc);
      if(delta_2 < a.back() && delta_2 >a.front())
	y2 = gsl_spline_eval(dist_table,delta_2,acc);
      double deriv = (y2 - y1)/(delta_2 - delta_1);
      double x_avg = (delta_2 + delta_1)/2;
      f->SetPoint(i,x_avg,deriv);
    }
  
  return f;
}


TGraph * GraphFunc(vector<double> &a, vector<double> &b){
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
void SetDataTable (const std::string& filename, vector<double> &x, vector<double> &y) {
  ifstream file (filename.c_str());
  double prev= 0;
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      x.push_back(d_energy/1000);
      y.push_back(d_value);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  return;
}

double calcInt(vector<double> &a, vector<double> &b, double delta){
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
  gsl_integration_qag (F_2, x1, delta, 0, 1e-2, 10000,3,w, &result, &error);    

  gsl_integration_workspace_free (w);
  gsl_interp_accel_free(acc);
  gsl_spline_free(dist_table);

  return result;
}

/*
TGraph * GraphCDF(vector<double> &a, vector<double> &b){
  int size_x = a.size();
  int npoints = 2*a.size();
  TGraph * f = new TGraph(npoints);
  //  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  //  gsl_spline *dist_table = gsl_spline_alloc(gsl_interp_akima,size_x);
  //  gsl_spline_init(dist_table,a.data(),b.data(),size_x);
  
  double x1 = a.front(), x2 = a.back();
  if(x1 < 1e-3) x1 = .1;//using log cant have nan errors
  double log_energy_step = (log(x2)-log(x1))/npoints;
  //  double scale = calcInt(a,b,a.back());

  for(int i = 0 ;i < npoints; i++)
    {
      double delta = i*log_energy_step + log(x1) ;
      delta = exp(delta);
      //      double value = 0;
      double value = calcInt(a,b,delta);
      f->SetPoint(i,delta,value);
      cout<<value<<endl;
    }

  return f;
}
      
*/

int main()
{
  Dielectric argon{1.66e-3,38.,18,"argon"};
  argon.SetPhotoCross("../argon_10_500_Sakamoto_ext.dat");

  //Set real and imaginary dielectric tables for faster computaiton
  if(argon.SetRealTable("real_argon.dat") == false)
    {
      argon.GetRealDielectric();
      argon.WriteToFile("real");
    }
  if(argon.SetImgTable("img_argon.dat") == false)
    {
      argon.GetImgDielectric();
      argon.WriteToFile("img");
    }  

  //  argon.GetRelCrossSection(3.16);
  //  TGraph *cross = argon.DrawCrossSection(true);
  //  TGraph *bichsel_0 = argon.DrawBichselSeg(938,1500,1.2,1e4,0,1e6,0);
  TGraph *bichsel_0 = argon.DrawBichselSeg(938,3470.6,1.2,1e4,0,1e6,0);
  TGraph *bichsel_1 = argon.DrawBichselSeg(938,3470.6,1.2,1e4,0,1e6,0);
  TH1D *c_dist = argon.GetMCdist(bichsel_0,50,.7,1e6);
  TH1D *c_dist_1 = argon.GetMCdist(bichsel_1,10,.7,1e6);

  TGraph *scaling = argon.GetCScaling(c_dist,c_dist_1);

  /*
    vector<double> cdf_x, cdf_y, dist_x,dist_y;
    SetDataTable("bichselStragCompare.dat",cdf_x,cdf_y);
    //    double scale = 2.19269;
    double scale = 2.25;
    for(int i = 0 ;i < cdf_x.size(); i++)
      cdf_y.at(i) /= scale;
    TGraph *cdf_calc = GraphCDF(cdf_x,cdf_y);
    TGraph *dist_paper = GraphFunc(cdf_x,cdf_y);
    
    SetDataTable("bichselStragCDF.dat",dist_x,dist_y);
    TGraph *cdf_from_paper = GraphFunc(dist_x,dist_y);//
    TGraph *dist_from_cdf = GraphCDF2Dist(dist_x,dist_y);
  */
  /*
  TCanvas *c7 = new TCanvas("c7","c7",1);
  c7 -> SetLogx();
  //  c7 -> SetLogy();
  //  bichsel -> GetYaxis() -> SetRangeUser(1e-9,2e-4);
  //  bichsel_0 -> GetXaxis() -> SetRangeUser(1e3,1e6);
  //  bichsel_0 -> GetXaxis() -> SetRangeUser(.1,1e3);
  //  bichsel_0 -> Draw();
  cdf_from_paper -> SetLineColor(4);
  cdf_from_paper -> Draw("");
  cdf_calc -> SetLineColor(2);
  cdf_calc -> Draw("same");


  c7 -> SaveAs("bichsel_seg.png");
  */
  TCanvas *c1 = new TCanvas("c1","c1",1);
  //  c1 -> SetLogx();
  //  c1 -> SetLogy();
  //  dist_from_cdf -> GetYaxis() -> SetRangeUser(0,.5);
  //  dist_from_cdf -> Draw();
  //  dist_paper->SetLineColor(2);
  //  dist_paper->Draw("same");
  //  bichsel_0->GetYaxis()->SetRangeUser(0,1);
  //  c_dist -> Draw("");
  //  c_dist -> GetXaxis() -> SetRangeUser(1e3,1e4);
  //  bichsel_0 -> SetLineColor(2);
  //  bichsel_0 -> Draw("same");
  //  scaling -> GetXaxis() -> SetRangeUser(2e3,3e3);
  scaling -> GetXaxis() -> SetRangeUser(1e3,2e3);
  scaling -> GetXaxis() -> SetRangeUser(1e3,2e3);
    scaling ->Draw("APO");

  c1 -> SaveAs("c_dist.png");
  
  return 0;
}
