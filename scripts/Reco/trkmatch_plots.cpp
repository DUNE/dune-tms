#include<string>
#include<stdio.h>
#include<cstdlib>
#include<iostream>
#include<math.h>
//#include<DUNEStyle.h>

using namespace std;

float delta_x(float XStart, float XEnd, float ZStart, float ZEnd, float XDir, float ZDir){
    float ExtrapX, XOffset, DeltaX;
    float dxdz;
	float ZGap = ZEnd - ZStart;

    dxdz = XDir / ZDir;

    XOffset = dxdz * ZGap;

    ExtrapX = XOffset + XStart;

    DeltaX = ExtrapX - XEnd;
   
    return DeltaX;
}

float delta_y(float YStart, float YEnd, float ZStart, float ZEnd, float YDir, float ZDir){
    float ExtrapY, YOffset, DeltaY;
    float dydz;
	float ZGap = ZEnd - ZStart;

    dydz = YDir / ZDir;

    YOffset = dydz *  ZGap;

	ExtrapY = YOffset + YStart;

    DeltaY = ExtrapY - YEnd;

    return DeltaY;
}

float delta_theta(float XDir1, float YDir1, float ZDir1, float XDir2, float YDir2, float ZDir2){
    float DeltaTheta, DirDot, abs1, abs2;

    DirDot = (XDir1 * XDir2) + (YDir1 * YDir2) + (ZDir1 * ZDir2);

    abs1 = sqrt((XDir1 * XDir1) + (YDir1 * YDir1) + (ZDir1 * ZDir1));
    abs2 = sqrt((XDir2 * XDir2) + (YDir2 * YDir2) + (ZDir2 * ZDir2));

    DeltaTheta = acos(DirDot / (abs1 * abs2)) * (180/3.14); 

    return DeltaTheta;

}

void trkmatch_plots(std::string filename){
	std::string filepath = "/pnfs/dune/persistent/users/kleykamp/tmsreco_combined_files/" + filename + ".tmsreco.root";
	std::unique_ptr<TFile> myFile(TFile::Open(filepath.c_str()));

	gStyle->SetOptStat(0);
        	  
	 if(!myFile || myFile->IsZombie()){
            std::cerr << "Error opening file" << endl;
            exit(-1);
	 }

	double mean, sigma;
    auto c1 = new TCanvas("");

	//declaring histograms
	TH1F *deltax_tms_truth = new TH1F("", "Delta X", 100, -600, 600);
    TH1F *deltay_tms_truth = new TH1F("", "Delta Y", 100, -600, 600);

	TH1F *deltax_tms_reco = new TH1F("", "Delta X, Reco", 100, -600, 600);
	TH1F *deltay_tms_reco = new TH1F("", "Delta Y, Reco", 100, -600, 600);

	TH1F *delta_theta_x_truth = new TH1F("", "", 80, 0, 35);
	TH1F *delta_theta_y_truth = new TH1F("", "", 80, 0, 35);

	TH1F *delta_theta_x_reco = new TH1F("", "", 80, 0, 35);
	TH1F *delta_theta_y_reco = new TH1F("", "", 80, 0, 35);


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
		truth->GetEntry(i);
		reco->GetEntry(i);

		int nTracksTrue = truth->GetLeaf("RecoTrackN")->GetValue(0);
		int nTracksReco = reco->GetLeaf("nTracks")->GetValue(0);

		int j, k, l;
		for(j=0, k=0, l=0; (l<= nTracksReco); j+=4, k+=3, l++){
			nhits_intms = reco->GetLeaf("nHits")->GetValue(l);
			pdg = truth->GetLeaf("PDG")->GetValue(l);

			//for some reason this logic works but just checking for != 13 || != -13 doesn't?
			if(pdg == 13 || pdg == -13){
			}
			else{
				//cout << "not muon, pdg " << pdg << endl;
				continue;
			}

			true_tms_start_x = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionEnteringTMS")->GetValue(j);
			true_tms_start_y = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionEnteringTMS")->GetValue(j+1);
			true_tms_start_z = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionEnteringTMS")->GetValue(j+2);
			
			true_lar_end_x = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionLeavingLAr")->GetValue(j);
			true_lar_end_y = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionLeavingLAr")->GetValue(j+1);
			true_lar_end_z = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionLeavingLAr")->GetValue(j+2);

			true_tms_p_start_x = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumEnteringTMS")->GetValue(j);
			true_tms_p_start_y = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumEnteringTMS")->GetValue(j+1);
			true_tms_p_start_z = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumEnteringTMS")->GetValue(j+2);

			true_lar_p_end_x = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumLeavingLAr")->GetValue(j);
			true_lar_p_end_y = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumLeavingLAr")->GetValue(j+1);
			true_lar_p_end_z = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumLeavingLAr")->GetValue(j+2);

			reco_tms_start_x = reco->GetLeaf("StartPos")->GetValue(k);
			reco_tms_start_y = reco->GetLeaf("StartPos")->GetValue(k+1);
			reco_tms_start_z = reco->GetLeaf("StartPos")->GetValue(k+2);

			reco_tms_startdir_x = reco->GetLeaf("StartDirection")->GetValue(k);
			reco_tms_startdir_y = reco->GetLeaf("StartDirection")->GetValue(k+1);
			reco_tms_startdir_z = reco->GetLeaf("StartDirection")->GetValue(k+2);

			float true_deathpos_x = truth->GetLeaf("DeathPosition")->GetValue(j);
			float true_deathpos_y = truth->GetLeaf("DeathPosition")->GetValue(j+1);
			float true_deathpos_z = truth->GetLeaf("DeathPosition")->GetValue(j+2);

			//double check that direction variables are normalized correctly
			float reco_dirlen = sqrt((reco_tms_startdir_x * reco_tms_startdir_x) + (reco_tms_startdir_y * reco_tms_startdir_y) + (reco_tms_startdir_z * reco_tms_startdir_z));
			reco_tms_startdir_x = reco_tms_startdir_x / reco_dirlen;
			reco_tms_startdir_y = reco_tms_startdir_y / reco_dirlen;
			reco_tms_startdir_z = reco_tms_startdir_z / reco_dirlen;

			//containment checks
			if(true_tms_start_x > 3300 || true_tms_start_x < -3300 || true_tms_start_y > 160 || true_tms_start_y < -2850){
				//cout << "did not end in tms" << endl;
				continue;
			}
			if(true_lar_end_x < -4500 || true_lar_end_x > 3700 || true_lar_end_y < -3200 || true_lar_end_y > 1000){
				//cout << "did not start in lar" << endl;
				continue;
			}
			if(true_deathpos_z > 18314){
				//cout << "punched out the back" << endl;
			    continue;
			}
			
			//cout << "passed cuts, yay!" << endl; //debugging
			
			//finding direction cosines from momentum information because true direction is not saved well
			float p_abs_tms_truth = sqrt((true_tms_p_start_x * true_tms_p_start_x) + (true_tms_p_start_y * true_tms_p_start_y) + (true_tms_p_start_z * true_tms_p_start_z));
			float dx_tms_truth = true_tms_p_start_x / p_abs_tms_truth;
			float dy_tms_truth = true_tms_p_start_y / p_abs_tms_truth;
			float dz_tms_truth = true_tms_p_start_z / p_abs_tms_truth;

			float p_abs_lar_truth = sqrt((true_lar_p_end_x * true_lar_p_end_x) + (true_lar_p_end_y * true_lar_p_end_y) + (true_lar_p_end_z * true_lar_p_end_z));
			float dx_lar_truth = true_lar_p_end_x / p_abs_lar_truth;
			float dy_lar_truth = true_lar_p_end_y / p_abs_lar_truth;
			float dz_lar_truth = true_lar_p_end_z / p_abs_lar_truth;

			
			float delta_theta_xz_tmstruth = delta_theta(dx_tms_truth, 0, dz_tms_truth, dx_lar_truth, 0, dz_lar_truth);
			float delta_theta_yz_tmstruth = delta_theta(0, dy_tms_truth, dz_tms_truth, 0, dy_lar_truth, dz_lar_truth);

			float delta_theta_xz_tmsreco = delta_theta(reco_tms_start_x, 0, reco_tms_start_z, dx_lar_truth, 0, dz_lar_truth);
			float delta_theta_yz_tmsreco = delta_theta(0, reco_tms_start_y, reco_tms_start_z, 0, dy_lar_truth, dz_lar_truth);

			//fill histograms
			deltax_tms_truth->Fill(delta_x(true_lar_end_x, true_tms_start_x, true_lar_end_z, true_tms_start_z, dx_lar_truth, dz_lar_truth));
			deltay_tms_truth->Fill(delta_y(true_lar_end_y, true_tms_start_y, true_lar_end_z, true_tms_start_z, dy_lar_truth, dz_lar_truth));

			deltax_tms_reco->Fill(delta_x(true_lar_end_x, reco_tms_start_x, true_lar_end_z, reco_tms_start_z, dx_lar_truth, dz_lar_truth));
			deltay_tms_reco->Fill(delta_y(true_lar_end_y, reco_tms_start_y, true_lar_end_z, reco_tms_start_z, dy_lar_truth, dz_lar_truth));

			delta_theta_x_truth->Fill(delta_theta_xz_tmstruth);
			delta_theta_y_truth->Fill(delta_theta_yz_tmstruth);

			delta_theta_x_reco->Fill(delta_theta_xz_tmsreco);
			delta_theta_y_reco->Fill(delta_theta_yz_tmsreco);
		}
	}

	//delta x and delta y plots
	deltax_tms_truth->SetTitle("\\Delta X (TMS Truth)");
	deltax_tms_truth->Draw();
	std::string printname = "plots/" + filename + "/deltax_tms_truth.png";
    c1->Print(printname.c_str());

	deltay_tms_truth->SetTitle("\\Delta Y (TMS Truth)");
	deltay_tms_truth->Draw();
	printname = "plots/" + filename + "/deltay_tms_truth.png";
	c1->Print(printname.c_str());

	deltax_tms_reco->SetTitle("\\Delta X (TMS Reco.)");
	deltax_tms_reco->Draw();
	printname = "plots/" + filename + "/deltax_tms_reco.png";
    c1->Print(printname.c_str());

	deltay_tms_reco->SetTitle("\\Delta Y (TMS Reco.)");
	deltay_tms_reco->Draw();
	printname = "plots/" + filename + "/deltay_tms_reco.png";
	c1->Print(printname.c_str());

	//delta theta_x and theta_y plots
	delta_theta_x_truth->SetTitle("\\Delta \\theta_x (TMS Truth)");
	delta_theta_x_truth->Draw();
	printname = "plots/" + filename + "/delta_theta_x_truth.png";
    c1->Print(printname.c_str());

	delta_theta_y_truth->SetTitle("\\Delta \\theta_y (TMS Truth)");
	delta_theta_y_truth->Draw();
	printname = "plots/" + filename + "/delta_theta_y_truth.png";
	c1->Print(printname.c_str());

	delta_theta_x_reco->SetTitle("\\Delta \\theta_x (TMS Reco.)");
	delta_theta_x_reco->Draw();
	printname = "plots/" + filename + "/delta_theta_x_reco.png";
    c1->Print(printname.c_str());

	delta_theta_y_reco->SetTitle("\\Delta \\theta_y (TMS Reco.)");
	delta_theta_y_reco->Draw();
	printname = "plots/" + filename + "/delta_theta_y_reco.png";
	c1->Print(printname.c_str());

	 
}
