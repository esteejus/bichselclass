void theory2graph(const std::string& filename){

  vector<double> theory_x, theory_y;
  ifstream file (filename.c_str());
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      theory_x.push_back(d_energy);
      theory_y.push_back(d_value);
    }
  }
  else std::cout << "Unable to open Straggling theory "<< filename << std::endl;


  TGraph *th = new TGraph(theory_x.size(),theory_x.data(),theory_y.data());
  th->Draw();

  return;

}
