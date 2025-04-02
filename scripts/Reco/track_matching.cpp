#include<DUNEStyle.h>
#include<iostream>
#include<string>
#include<cstdlib>

using namespace std;

float ZGap = 2250; //the distance between ND-LAr and TMS in mm

bool isTMSContained(float xpos, float ypos, float zpos){
	bool out = true;
	
	if(xpos > 3300 || xpos < -3300){
		out = false;
	}
	if(ypos > 160 || ypos < -2850){
		out = false;
	}
	if(zpos > 18314 || zpos < 11362){ //tms end is 18314mm for new thick/thin geometry; is 13500mm for thin planes only
		out = false;
	}

	return out;
}

bool isMuon(int PDG){
	bool out = false;

	if(std::abs(PDG) == 13){
		out = true;
	}

	return out;
}

int main(std::string filepath){

	std::unique_ptr<TFile> myFile(TFile::Open(filepath.c_str()));
	
	if(!myFile || myFile->IsZombie()){
		std::cerr << "Error opening file" << endl;
		exit(-1);
	}

	TTree *truth = myFile->Get<TTree>("Truth_Info");
	TTree *reco = myFile->Get<TTree>("Reco_Tree");

	int nentries = truth->GetEntries();

	return 0;
}
