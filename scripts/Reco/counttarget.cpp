#include <algorithm>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"

void counttarget() {
  TFile *f = new TFile("neutrino_merge_0p1Ecut.root");
  TTree *t = (TTree*)f->Get("Truth_Info");

  std::string *Interaction;
  float Muon_Vertex[4];
  t->SetBranchStatus("*", false);
  t->SetBranchStatus("Interaction", true);
  t->SetBranchAddress("Interaction", &Interaction);
  t->SetBranchStatus("Muon_Vertex", true);
  t->SetBranchAddress("Muon_Vertex", Muon_Vertex);

  // Just get the interaction string
  int nIron = 0;
  int nCarbon = 0;
  int nHydrogen = 0;
  int nNitrogen = 0;
  int nSilicon = 0;
  int nOthers = 0;

  int nEntries = t->GetEntries();
  //nEntries = 20000;
  std::vector<std::string> targets;
  for (int i = 0; i < nEntries; ++i) {
    t->GetEntry(i);
    if (std::find(targets.begin(), targets.end(), *Interaction) != targets.end()) targets.push_back(*Interaction);
    if ((*Interaction).find("tgt:1000260560") != std::string::npos) nIron++;
    else if ((*Interaction).find("tgt:1000060120") != std::string::npos) nCarbon++;
    else if ((*Interaction).find("tgt:1000140280") != std::string::npos) nSilicon++;
    else if ((*Interaction).find("tgt:1000070140") != std::string::npos) nNitrogen++;
    else if ((*Interaction).find("tgt:1000010010") != std::string::npos) nHydrogen++;
    else {
      std::cout << *Interaction << std::endl;
      nOthers++;
    }

  }
  std::cout << "Iron: " << nIron << "/" << nEntries << " (" << double(nIron)/nEntries*100. << "%)" << std::endl;
  std::cout << "Carbon: " << nCarbon << "/" << nEntries << " (" << double(nCarbon)/nEntries*100. << "%)" << std::endl;
  std::cout << "Nitrogen: " << nNitrogen << "/" << nEntries << " (" << double(nNitrogen)/nEntries*100. << "%)" << std::endl;
  std::cout << "Silicon: " << nSilicon << "/" << nEntries << " (" << double(nSilicon)/nEntries*100. << "%)" << std::endl;
  std::cout << "Hydrogen: " << nHydrogen << "/" << nEntries << " (" << double(nHydrogen)/nEntries*100. << "%)" << std::endl;
  std::cout << "Other: " << nOthers << "/" << nEntries << " (" << double(nOthers)/nEntries*100. << "%)" << std::endl;

  std::cout << "numebr of targets: " << targets.size() << std::endl;

}
