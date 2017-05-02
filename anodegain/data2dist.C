void setdistfile(const std::string& filename,TH1D *data) {
  ifstream file (filename.c_str());
  double d_dedx=-99.;
  //  string rnd_num;
  std::string index,line;

  if (file.is_open()){
    getline(file,line); //skip header
    while(getline (file,line) ){
      
      std::istringstream in(line);
      in>>d_dedx;
      data->Fill(d_dedx);
    }
  }
  else std::cout << "Unable to open Inverted Cross seciton file with SetInvXSec" << std::endl;

  data->Scale(1./data->Integral("width"));
  
  return;
}

void data2dist(string name){
  TH1D *dist = new TH1D("dist","dist",2000,0,2000);
  setdistfile(name,dist);
  ofstream out;
  out.open("anode_dist.dat");

  //out<<"Particle "<<name<<" momenum [MeV/c] "<< mom << " mass [MeV/c^2] "<<mass;
  int num_binsx = dist->GetNbinsX();
  for(int i = 1; i<= num_binsx; i++){
    out<<dist->GetBinCenter(i)<<"\t"<<dist->GetBinContent(i)<<endl;  
  }
  out.close();  
 




  return;
}
