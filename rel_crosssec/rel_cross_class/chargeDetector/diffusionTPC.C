


void diffusionTPC()
{


  double eVtoelectron = 26.4;// eV/electron (actually ion pair)
  double segment_l = 1.//cm
  double delatE = 1.2; //keV/cm




  TRandom3 *ran = new TRandom3(0);
  double diff_coeff = .024; //sqrt(cm)-1 transverse diffusion
  double drift_l = 40;// drift length [cm]
  double sigmaT = diff_coeff * sqrt(drift_l);

  int num_electrons = down(deltaE/eVtoelctron);//

  vector<double> x,y,z,charge;
  
  for(int iEl = 0 ;iEl < num_electrons; iEl ++)
    {


    }
  
  double dr = ran -> Gaus(0,sigmaT);
  double angle = ran -> Uniform(2*TMath::Pi());
  double dx = dr * cos(angle);
  double dz = dr * sin(angle);

  TF1 *diff_fcn = new TF1("diff_fcn","[0]*TMath::Gaus(x,0,[1])",-1,1);
  cout<<sigmaT<<endl;
  diff_fcn->SetParameters(1,sigmaT);
  double scale = diff_fcn->Integral(-100,100);
  diff_fcn->SetParameter(0,1./scale);

  diff_fcn->Draw();


  return;
}
