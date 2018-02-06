/*#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>
#include "TF1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "../dielectric.h"
#include "TRandom3.h"

using namespace std;
*/

void plotTable()//int main()
{
  TFile *f = new TFile("analyzedtable.root");
  TTree *t = (TTree *)f->Get("tree");

  double bg, slope, offset;
  TH1D *refdist = nullptr;
  TGraph *scaledref = nullptr;
  TGraph *scalingrel = nullptr;
  TH1D *cdist = nullptr;
  
  t->SetBranchAddress("bg",&bg);
  t->SetBranchAddress("slope",&slope);
  t->SetBranchAddress("offset",&offset);
  t->SetBranchAddress("refdist",&refdist);
  t->SetBranchAddress("scaledref",&scaledref);
  t->SetBranchAddress("scalingrel",&scalingrel);
  t->SetBranchAddress("cdist",&cdist);

  for(int iEntry = 0;iEntry < t->GetEntries(); iEntry++)
    {
      t->GetEntry(iEntry);
      cout<<bg<<" "<<slope<<" "<<offset<<endl;
      /*
      TCanvas *c1 = new TCanvas("c1","c1",1);
      c1->SetLogy();
      c1->SetLogx();
      cdist->Draw();
      scaledref->SetMarkerStyle(20);
      scaledref->SetMarkerSize(.5);
      scaledref->Draw("same PO");
      c1->SaveAs(Form("%i.png",iEntry));
      */
    }

  return 0;
}
