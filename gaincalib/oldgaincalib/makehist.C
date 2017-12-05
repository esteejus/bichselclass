

void makehist()
{
  TH1D *dist = new TH1D("olddist","olddist",1000,0,1000);
 ifstream file ("data_dedx_p.dat");
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      dist->Fill(d_energy);
    }
  }
    TFile *out = new TFile("olddata.root","RECREATE");
    dist->Write();


  return;
}
