#include<string>
#include<stdio.h>
#include<cstdlib>
#include<iostream>
#include<math.h>
#include<DUNEStyle.h>

using namespace std;

float ZGap = 2250;
float tmstolar = -1;
float lartotms = 1;

float delta_x(float XStart, float XEnd, float XDir, float ZDir, float Direction){
    float ExtrapX, XOffset, DeltaX;
    float dxdz;

    dxdz = XDir / ZDir;

    XOffset = dxdz * ZGap * Direction;

    ExtrapX = XOffset + XStart;

    DeltaX = ExtrapX - XEnd;
   
    return DeltaX;
}

float delta_y(float YStart, float YEnd, float YDir, float ZDir, float Direction){
    float ExtrapY, YOffset, DeltaY;
    float dydz;

    dydz = YDir / ZDir;

    YOffset = dydz * Direction * ZGap;

	ExtrapY = YOffset + YStart;

    DeltaY = ExtrapY - YEnd;

    return DeltaY;
}

float delta_r(float XStart, float YStart, float XEnd, float YEnd, float XDir, float YDir, float ZDir, float Direction){
    float DeltaX, DeltaY, DeltaR;

    DeltaX = delta_x(XStart, XEnd, XDir, ZDir, Direction);
    DeltaY = delta_y(YStart, YEnd, YDir, ZDir, Direction);

    DeltaR = sqrt(pow(DeltaX, 2) + pow(DeltaY, 2));

    return DeltaR;
}

float delta_theta(float XDir1, float YDir1, float ZDir1, float XDir2, float YDir2, float ZDir2){
    float DeltaTheta, DirDot, abs1, abs2;

    DirDot = (XDir1 * XDir2) + (YDir1 * YDir2) + (ZDir1 * ZDir2);

    abs1 = sqrt((XDir1 * XDir1) + (YDir1 * YDir1) + (ZDir1 * ZDir1));
    abs2 = sqrt((XDir2 * XDir2) + (YDir2 * YDir2) + (ZDir2 * ZDir2));

    DeltaTheta = acos(DirDot / (abs1 * abs2)) * (180/3.14); 

    return DeltaTheta;

}

void track_matching_clean(std::string filename){
	std::unique_ptr<TFile> myFile(TFile::Open(filename.c_str()));

	gStyle->SetOptStat(0);
        	  
	 if(!myFile || myFile->IsZombie()){
            std::cerr << "Error opening file" << endl;
            exit(-1);
	 }

	 //defining truth variables
	 float true_tms_start_x, true_tms_start_y, true_tms_start_z, true_lar_end_x, true_lar_end_y, true_lar_end_z;
	 float true_tms_p_start_x, true_tms_p_start_y, true_tms_p_start_z, true_lar_p_end_x, true_lar_p_end_y, true_lar_p_end_z;
	 int pdg;

	 //defining reco variables
	 float reco_tms_start_x, reco_tms_start_y, reco_tms_start_z, reco_lar_end_x, reco_lar_end_y, reco_lar_end_z;
	 float reco_tms_startdir_x, reco_tms_startdir_y, reco_tms_startdir_z, reco_lar_enddir_x, reco_lar_enddir_y, reco_lar_enddir_z;
	 float nhits_intms;

	TTree *truth = myFile->Get<TTree>("Truth_Info");
	TTree *reco = myFile->Get<TTree>("Reco_Tree");

	int nentries = truth->GetEntries();

	for(int i = 0; i <= nentries; i++){
		truth-GetEntry(i);
		reco->GetEntry(i);

		int nTracksTrue = truth->GetLeaf("RecoTrackN")->GetValue(0);
		int nTracksReco = reco->GetLeaf("nTracks")->GetValue(0);

		
	}

	 
}
