#include <iostream>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <vector>
#include "TCanvas.h"
#include "TF2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1.h"

using namespace std;


int main()
{
  //estimate 2D Guass
  double xmin = -20;
  double xmax = 20;
  double ymin = -20;
  double ymax = 20;

  TF2 *fcn = new TF2("gaus","TMath::Gaus(x,0,1,kTRUE)*TMath::Gaus(y,0,2.5,kTRUE)",xmin,xmax,ymin,ymax);
  //  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  int xgrid_points = 1e2;
  int ygrid_points = 1e2;
  double xstep = (xmax - xmin)/(xgrid_points-1);
  double ystep = (ymax - ymin)/(ygrid_points-1);

  double xa[xgrid_points];
  double ya[ygrid_points];
  double za[xgrid_points*ygrid_points];

  size_t nx = sizeof(xa) / sizeof(xa[0]);
  size_t ny = sizeof(ya) / sizeof(ya[0]);

  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  for(int i = 0; i < xgrid_points; i++)
    {
      for(int j = 0; j < ygrid_points; j++)
	{
	  double x =  xstep * i + xmin;
	  double y =  ystep * j + ymin;
	  double z = fcn->Eval(x,y);
	  xa[i] = x;
	  ya[j] = y; 
	  gsl_spline2d_set(spline, za, i, j, z);
	  cout<<x<<" "<<y<<endl;
	}
    }
  
  gsl_spline2d_init(spline, xa, ya, za, nx, ny);

  TRandom3 *ran = new TRandom3(0);
  TH1D *residual =  new TH1D("residual","residual",1000,-1,1);
  for(int i = 0; i < 1e4; i++)
    {
      double x = ran -> Uniform(xmin,xmax); 
      double y =  ran -> Uniform(ymin,ymax);
      double z = fcn->Eval(x,y);
      double est = gsl_spline2d_eval(spline,x,y,xacc, yacc);
      //      cout<<"difference "<<est-z<<endl;
      residual->Fill((est-z)/z);
    }

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  TCanvas *c1 = new TCanvas("c1","c1",1);
  residual->Draw();
  c1->SaveAs("bicubic_res.png");
  
  return 0;

}
