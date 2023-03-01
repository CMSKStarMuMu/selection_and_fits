#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>       


void add_composed_weight(int year = 2016, int type = 0, int mc = 0){
  
  // define input files based on the year and the dataset
  std::string version = "XGBv8";

  std::string     channel = "LMNR";
  if (type == 1)  channel = "JPSI"; 
  if (type == 2)  channel = "PSI"; 

  std::string     f_channel = "LMNR";
  if (type == 1)  f_channel = "Jpsi"; 
  if (type == 2)  f_channel = "Psi"; 

  std::string       string_nonan = "";
  if (year == 2017) string_nonan = "_noNan"; 

  std::cout << "creating new tree for year: " << year << "  channel: " << channel << std::endl;

  auto f_input = TFile::Open(Form("/eos/cms/store/group/phys_bphys/fiorendi/p5prime/ntuples/after_nominal_selection/%dMC_%s_noIP2D%s_addxcutvariable.root", year, channel.c_str(), string_nonan.c_str() ));
  auto tree = (TTree*)f_input->Get("ntuple");

  std::cout << "input file: " << f_input->GetName() << std::endl;
  std::cout << "input entries: " << tree->GetEntries() << std::endl;

  std::string friend_name = Form("/eos/user/a/aboletti/BdToKstarMuMu/fileIndex/MC-%s-%s/%d.root", f_channel.c_str(), version.c_str(), year);
  tree->AddFriend("wTree", friend_name.c_str());
  std::cout << "friend file: " << friend_name.c_str() << std::endl;

  TFile* out_file = new TFile(Form("input_with_weights/%dMC_%s_noIP2D%s_addxcutvariable_withMCw_%s.root", year, channel.c_str(), string_nonan.c_str(), version.c_str()),"RECREATE");
  std::cout << "output file : " << out_file -> GetName() << std::endl;
  
  float MCw, weight;
  float tot_weight;

  tree->SetBranchStatus("*",1);
  tree->SetBranchAddress("MCw", &MCw);
  tree->SetBranchAddress("weight", &weight);

  int nentries = tree->GetEntries();
  std::cout << "original n entries: " << nentries << std::endl;

  out_file -> cd();
  TTree* out_tree = (TTree*) tree->CloneTree(0);
  out_tree->Branch("tot_weight", &tot_weight);
  out_tree->SetAutoSave(1000000);
 
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent = tree-> GetEntry(eventNo);
    tot_weight = MCw * weight;
    out_tree -> Fill();
  }

  std::cout << "final n entries: " << out_tree->GetEntries() << std::endl;
  out_file -> cd();
  out_tree -> Write();

  return;
}