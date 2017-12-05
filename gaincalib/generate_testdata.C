



void generate_testdata()
{
  double slope = 15.34;
  double offset = 10.33;

  double mean = 0.;
  double sigma = .1;
  double norm = (1./sqrt(2*TMath::Pi()*sigma*sigma));
  int events = 1e4;

  TF1 * theory = new TF1("theory","[0]*TMath::Gaus(x,[1],[2])",-20,20);
  theory -> SetParameters(norm,mean,sigma);

  TRandom3 *ran = new TRandom3(0);



  TH1D *data = new TH1D("data","data",200,0,20);
  for(int i = 0;i < events ;i++)
    {
      double exp_value = ran->Gaus(mean,sigma);
      exp_value = exp_value*slope + offset;
      data->Fill(exp_value);
    }
  data->Scale(1./data->Integral("width"));
  
  cout<<"norm is "<<  data->Integral("width")<<endl;
  ofstream outdata;
  ofstream outtheory;

  outdata.open("test_data.dat");
    for(int i = 1;i < data->GetXaxis()->GetNbins()+1;i++)
      {
	outdata<<data->GetBinCenter(i)<<"\t"<<data->GetBinContent(i)<<endl;
      }

    int nstep = 1e4;
    double step_size = 40./nstep;
    outtheory.open("test_theory.dat");
    for(int i = 0; i<nstep;i++)
      {
	double x = i*step_size - 20.;
	double value = theory->Eval(x);
	if(value > 1e-10)
	  outtheory<<x<< "\t"<<value<<endl;
      }

    data->Draw();
    
  return;
}
