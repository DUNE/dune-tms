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

bool isLArContained(float xpos, float ypos, float zpos){
	bool out = true;

	if(xpos > 3700 || xpos < -4500){
		out = false;
	}
	if(ypos > 1000 || ypos < -3200){
		out = false;
	}
	if(zpos > 9200 || zpos < 4100){
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

bool hasMinHits(int nHits){
	bool out = false;

	if(nHits >= 5){
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

	for(int i = 0; i <= nentries; i++){
		truth->GetEntry(i);
		reco->GetEntry(i);
		
		int nTracksTrue = truth->GetLeaf("RecoTrackN")->GetValue(0);
		
		for(int itrack = 0; itrack <= nTracksTrue; itrack++){
			auto end_array = truth.RecoTrackPrimaryParticleTruePositionTrackEnd[itrack];
			auto start_array = truth.RecoTrackPrimaryParticleTruePositionTrackStart[itrack];			

			float XTrackEndpoint = end_array[0];
			float YTrackEndpoint = end_array[1];
			float ZTrackEndpoint = end_array[2]; 
			
			float XTrackStartpoint = start_array[0];
			float YTrackStartPoint = start_array[1];
			float ZTrackStartPoint = start_array[2];		
	
			bool is_muon = isMuon(truth->GetLeaf("PDG")->GetValue(itrack));
			bool is_TMS_contained = isTMSContained(XTrackEndpoint, YTrackEndpoint, ZTrackEndpoint);		
			bool is_LAr_contained = isLArContained(XTrackStartpoint, YTrackStartpoint, ZTrackStartpoint);
			bool has_min_hits = hasMinHits(reco->GetLeaf("nHits")->GetValue(itrack));	

			if(is_muon && is_TMS_contained && is_LAr_contained){

			}
		}
	}

	return 0;
}
