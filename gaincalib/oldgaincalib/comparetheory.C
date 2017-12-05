

void comparetheory()
{

  vector<double> old_x, old_y;
  vector<double> new_x, new_y;
 ifstream file ("p_dedx_compare.data");
  double d_value=0,d_energy=0;
  std::string index,line;

  if (file.is_open()){
    while(getline (file,line) ){
      std::istringstream in(line);
      in>>d_energy;
      in>>d_value;
      old_x.push_back(d_energy);
      old_y.push_back(d_value);
    }
  }

  //ifstream file2 ("../bichsel_full_p_903_108.data");
  ifstream file2 ("./bichsel_p10_theory.dat");
 d_value=0,d_energy=0;
  std::string index2,line2;

  if (file2.is_open()){
    while(getline (file2,line2) ){
      std::istringstream in(line2);
      in>>d_energy;
      in>>d_value;
      
      new_x.push_back(d_energy);
      new_y.push_back(d_value);
    }
  }

  TGraph *new_g = new TGraph(new_x.size(),new_x.data(),new_y.data());
  TGraph *old_g = new TGraph(old_x.size(),old_x.data(),old_y.data());


  double new_scale = new_g->Integral(1,new_g->GetN());
  for(int i = 0 ; i < new_y.size();i++)
    new_y.at(i) = new_y.at(i)/new_scale;
  
  TGraph *scale_new = new TGraph(new_x.size(),new_x.data(),new_y.data());

  double scale = old_g->Integral(1,old_g->GetN());
  for(int i = 0 ; i < old_y.size();i++)
    old_y.at(i) = old_y.at(i)/scale;
  
  TGraph *scale_old = new TGraph(old_x.size(),old_x.data(),old_y.data());

  cout<<"Integral is now for new "<<scale_new->Integral(1,scale_new->GetN())<<endl;
  cout<<"Integral is now for old "<<scale_old->Integral(1,scale_old->GetN())<<endl;
  TCanvas *c1 = new TCanvas("c1","c1",1);
  c1->SetLogx();
  //  cout<<  new_g->Integral(1,new_g->GetN())<<endl;
  scale_old->SetLineColor(2);
  scale_new->Draw("ALO");
  scale_old->Draw("same LO");

  return;
}
