#include<DUNEStyle.h>
#include<iostream>
#include<string>
#include<cstdlib>

using namespace std;

float ZGap = 2250; //the distance between ND-LAr and TMS in mm

int main(std::string filepath){

	std::unique_ptr<TFile> myFile(TFile::Open(filepath.c_str()));
	
	TTree *truth = myFile->Get<TTree>("Truth_Info");
	TTree *reco = myFile->Get<TTree>("Reco_Tree");

	int nentries = truth->GetEntries();

	return 0;
}
