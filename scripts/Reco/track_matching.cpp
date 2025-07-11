//script to do basic track matching between tms and nd lar

#include<DUNEStyle.h>
#include<TFile.h>
#include<TH1D.h>

#include<iostream>
#include<string>
#include<cstdlib>
#include<stdio.h>
#include<math.h>

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

float delta_1d(float startpos, float endpos, float dir, float zdir, float trackdir){
	float offset, extrapolation, delta;	
	
	offset = (dir / zdir) * ZGap * trackdir;
	extrapolation = offset + startpos; 

	delta = extrapolation - endpos; 

	return delta;
}

float delta_theta(float xdir1, float ydir1, float zdir1, float xdir2, float ydir2, float zdir2){
	float direction_dot, abs_vec1, abs_vec2, delta_theta;
	
	direction_dot = (xdir1 * xdir2) + (ydir1 * ydir2) + (zdir1 * zdir2);
	abs_vec1 = sqrt((xdir1 * xdir1) + (ydir1 * ydir1) + (zdir1 * zdir1));
	abs_vec2 = sqrt((xdir2 * xdir2) + (ydir2 * ydir2) + (zdir2 * zdir2));

	delta_theta = acos(direction_dot / (abs_vec1 * abs_vec2)) * (180/3.14);
	
	return delta_theta;
}

bool ismatched(){
	bool is_matched = false;



	return is_matched;
}

int main(std::string inputFilename){

	std::unique_ptr<TFile> myFile(TFile::Open(inputFilename.c_str()));
	
	if(!myFile || myFile->IsZombie()){
		std::cerr << "Error opening file" << endl;
		exit(-1);
	}

	std::string directoryPath ="/exp/dune/data/users/"+ std::string(getenv("USER")) + "/dune-tms/Reco/Track_Matching_Plots/";

   	if(createDirectory(directoryPath)) {
       		std::cout << "Directory created: " << directoryPath << std::endl;
    	}
	else{
        	std::cerr << "Failed to create directory" << std::endl;
    	}

	// Create output filename
	std::string outputFilename = directoryPath + getOutputFilename(inputFilename);
		
	// Create TFile with the output filename
	TFile outputFile(outputFilename.c_str(), "RECREATE");	
		
	TTree *truth = myFile->Get<TTree>("Truth_Info");
	TTree *reco = myFile->Get<TTree>("Reco_Tree");

	float TrueTMSStart;
	float PDG;

	truth->SetBranchAddress("RecoTrackPrimaryParticleTruePositionEnteringTMS", &TrueTMSStart);
	truth->SetBranchAddress("PDG", PDG);

	int nentries = truth->GetEntries();

	for(int i = 0; i <= nentries; i++){
		truth->GetEntry(i);
		reco->GetEntry(i);
		
		int nTracksTrue = truth->GetLeaf("RecoTrackN")->GetValue(0);
			
		for(int itrack = 0; itrack <= nTracksTrue; itrack++){
			bool is_muon = isMuon(truth->GetLeaf("PDG")->GetValue(itrack));
			bool is_TMS_contained = isTMSContained(XTrackEndpoint, YTrackEndpoint, ZTrackEndpoint);		
			bool is_LAr_contained = isLArContained(XTrackStartpoint, YTrackStartpoint, ZTrackStartpoint);
			bool has_min_hits = hasMinHits(reco->GetLeaf("nHits")->GetValue(itrack));	

			if(is_muon && is_TMS_contained && is_LAr_contained && has_min_hits){
				cout << XTrackEndpoint << " X" << endl;
			}
		}
	}
	
	outputFile.Close();
	return 0;
}
