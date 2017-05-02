

void combinepid(){

  TFile *f = TFile::Open("pid_high_test.root");
  TH2D *high = (TH2D*)f->Get("pid");
  high->Draw("colz");

  double mass = 140.;

  int nbins_x = high->GetXaxis()->GetNbins();
  int nbins_y = high->GetYaxis()->GetNbins();

  TH2D *inv = new TH2D("inv","inv",nbins_y,0,16,nbins_x,

  for(int i = 0;i < nbins_x;++i){
    for(int j = 0;j< nbins_y;++j){
   
      int counts = high -> GetBinContent(i,j);
      double dedx = high->GetYaxis()->GetBinCenter(i);
      double mom = high->GetXaxis()->GetBinCenter(i);
 


  }
  }

  return;
}
