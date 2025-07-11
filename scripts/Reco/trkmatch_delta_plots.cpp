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
    std::string filepath = "root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr/dune/persistent/users/kleykamp/tmsreco_combined_files/" + filename + ".tmsreco.root";

    std::unique_ptr<TFile> myFile(TFile::Open(filepath.c_str()));

	gStyle->SetOptStat(0);
        	  
	 if(!myFile || myFile->IsZombie()){
            std::cerr << "Error opening file" << endl;
            exit(-1);
    }

    //defining input histograms
    /*
    TH1F *x_pos_tms_truth_h = new TH1F("","RecoTrackPrimaryParticleTruePositionEnteringTMS, X", 100, -3200, -3200);
    TH1F *y_pos_tms_truth_h = new TH1F("","RecoTrackPrimaryParticleTruePositionEnteringTMS, Y", 100, -3200, -3200);
    TH1F *z_pos_tms_truth_h = new TH1F("","RecoTrackPrimaryParticleTruePositionEnteringTMS, Z", 100, -3200, -3200);

    TH1F *x_pos_lar_truth_h = new TH1F("","RecoTrackPrimaryParticleTruePositionLeavingLAr, X", 100, -3200, -3200);
    TH1F *y_pos_lar_truth_h = new TH1F("","RecoTrackPrimaryParticleTruePositionLeavingLAr, Y", 100, -3200, -3200);
    TH1F *z_pos_lar_truth_h = new TH1F("","RecoTrackPrimaryParticleTruePositionLeavingLAr, Z", 100, -3200, -3200);

    TH1F *dx_tms_true_h = new TH1F("X Momentum at TMS Start / |P_TMS|","", 100, -1, 1);
    TH1F *dy_tms_true_h = new TH1F("Y Momentum at TMS Start / |P_TMS|","", 100, -1, 1);
    TH1F *dz_tms_true_h = new TH1F("Z Momentum at TMS Start / |P_TMS|","", 100, -1, 1);

    TH1F *dx_lar_true_h = new TH1F("X Momentum at LAr End / |P_LAr|","", 100, -1, 1);
    TH1F *dy_lar_true_h = new TH1F("Y Momentum at LAr End / |P_LAr|","", 100, -1, 1);
    TH1F *dz_lar_true_h = new TH1F("Z Momentum at LAr End / |P_LAr|","", 100, -1, 1); 

    TH1F *x_pos_tms_reco_h = new TH1F("", "TMS Start (X) Pos, Reco", 100, -3200, -3200);
    TH1F *y_pos_tms_reco_h = new TH1F("", "TMS Start (Y) Pos, Reco", 100, -3200, -3200);
    TH1F *z_pos_tms_reco_h = new TH1F("", "TMS Start (Z) Pos, Reco", 100, -3200, -3200);

    /* currently unneeded
    TH1F *x_pos_lar_reco_h = new TH1F("", "LAr End (X) Pos, Reco", 100, -3200, -3200);
    TH1F *y_pos_lar_reco_h = new TH1F("", "LAr End (Y) Pos, Reco", 100, -3200, -3200);
    TH1F *z_pos_lar_reco_h = new TH1F("", "LAr End (Z) Pos, Reco", 100, -3200, -3200);

    */
    /*
    TH1F *x_dir_tms_reco_h = new TH1F("", "TMS Start (X) Dir, Reco", 100, -1, 1);
    TH1F *y_dir_tms_reco_h = new TH1F("", "TMS Start (Y) Dir, Reco", 100, -1, 1);
    TH1F *z_dir_tms_reco_h = new TH1F("", "TMS Start (Z) Dir, Reco", 100, -1, 1);
    
    //defining delta histos
    TH1F *xhisto_tmstolar_true = new TH1F("", "DeltaX", 200, -1000, 1000);
    TH1F *yhisto_tmstolar_true = new TH1F("", "DeltaY", 200, -1000, 1000);
    TH1F *rhisto_tmstolar_true = new TH1F("", "DeltaR", 200, 0, 2000);
    
    TH1F *xhisto_lartotms_true = new TH1F("", "", 200, -1000, 1000);
    TH1F *yhisto_lartotms_true = new TH1F("", "", 200, -1000, 1000);
    TH1F *rhisto_lartotms_true = new TH1F("", "", 200, 0, 2000);

    TH1F *xhisto_lartotms_truetoreco = new TH1F("", "", 200, -1000, 1000);
    TH1F *yhisto_lartotms_truetoreco = new TH1F("", "", 200, -1000, 1000);
    TH1F *rhisto_lartotms_truetoreco = new TH1F("", "", 200, 0, 2000);

    TH1F *xhisto_tmstolar_recototrue = new TH1F("", "", 200, -1000, 1000);
    TH1F *yhisto_tmstolar_recototrue = new TH1F("", "", 200, -1000, 1000);
    TH1F *rhisto_tmstolar_recototrue = new TH1F("", "", 200, 0, 2000);

    TH1F *thetahisto_tms_truevreco = new TH1F("", "", 180, 0, 30);
    TH1F *thetahisto_tmslar_truth = new TH1F("", "", 180, 0, 30);
    TH1F *thetahisto_tmslar_recotrue = new TH1F("", "", 180, 0, 30);

    TH2F *deltay_vs_tracklength = new TH2F("", "", 33, -3200, -3200, 33, -3200, -3200);
    TH2F *y_dir_tms_reco_vs_tracklength = new TH2F("", "", 33, -3200, -3200, 33, -3200, -3200);
    TH2F *delta_theta_tmslar_recotrue_vs_tracklength = new TH2F("", "", 30, 0, 50, 30, 0, 7000);

    TH2F *theta_xz_tmsreco_vs_theta_xz_lartruth = new TH2F("", "", 60, -30, 30, 60, -30, 30);
    TH2F *theta_yz_tmsreco_vs_theta_yz_lartruth = new TH2F("", "", 60, -30, 30, 60, -30, 30);

    TH1F *delta_xz_tms_truereco = new TH1F("", "", 80, -25, 25);
    TH1F *delta_yz_tms_truereco = new TH1F("", "", 80, -25, 25);

    TH1F *delta_xz_tmslar_true = new TH1F("", "", 80, -10, 10);
    TH1F *delta_yz_tmslar_true = new TH1F("", "", 80, -10, 10);

    TH1F *delta_xz_tmslar_recotrue = new TH1F("", "", 80, -25, 25);
    TH1F *delta_yz_tmslar_recotrue = new TH1F("", "", 80, -25, 25);
    */

    TH1F *energyhisto = new TH1F("", "", 100, -1000, -1000);
    TPaveText *fitinfo = new TPaveText(.6, .65, .8, .85, "NB NDC");

    double mean, sigma;

    auto c1 = new TCanvas("");

    //defining truth variables
    Float_t XPosTMSStartTrue, YPosTMSStartTrue, ZPosTMSStartTrue, XPosLArEndTrue, YPosLArEndTrue, ZPosLArEndTrue;
    Float_t XMomTMSStartTrue, YMomTMSStartTrue, ZMomTMSStartTrue, XMomLArEndTrue, YMomLArEndTrue, ZMomLArEndTrue;
    Float_t XEndTrue, YEndTrue, ZEndTrue, XStartTrue, YStartTrue, ZStartTrue;

    //defining reco variables
    float XPosTMSStartReco, YPosTMSStartReco, ZPosTMSStartReco;
    float XDirTMSStartReco, YDirTMSStartReco, ZDirTMSStartReco;
    float TrackLengthInTMS, nHits;
    int PDG;

    TTree *truth = myFile->Get<TTree>("Truth_Info");
    TTree *reco = myFile->Get<TTree>("Reco_Tree");

    int nentries = truth->GetEntries();	

    for(int i = 0; i <= nentries; i++){
        truth->GetEntry(i);
        reco->GetEntry(i);

        int nTracksTrue = truth->GetLeaf("RecoTrackN")->GetValue(0);
        int nTracksReco = reco->GetLeaf("nTracks")->GetValue(0);

        if(nTracksTrue != nTracksReco){
            cout << "The number of tracks in Reco_Tree and Truth_Info does not match. nTracksReco = " << nTracksReco << " and nTracksTruth = " << nTracksTrue << endl;
            cout << "Event #: " << reco->GetLeaf("EventNo")->GetValue(0) << endl;
            continue; 
        }

        if(nTracksTrue == 0){
            continue;
        }

        int j, k, l;

        for(j = 0, k = 0, l = 0; (l <= nTracksReco); j+=4, k+=3, l++){
            TrackLengthInTMS = truth->GetLeaf("RecoTrackPrimaryParticleTrueTrackLengthInTMS")->GetValue(l);
            nHits = reco->GetLeaf("nHits")->GetValue(l);
            PDG = truth->GetLeaf("PDG")->GetValue(l);

            if(PDG == 13 || PDG == -13){
                //cout << "PDG: " << PDG << endl;
            }
            else{
                continue;
            }

            //pulling truth variables
            XPosTMSStartTrue = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionEnteringTMS")->GetValue(j);
            YPosTMSStartTrue = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionEnteringTMS")->GetValue(j+1);
            ZPosTMSStartTrue = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionEnteringTMS")->GetValue(j+2);

            XPosLArEndTrue = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionLeavingLAr")->GetValue(j);
            YPosLArEndTrue = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionLeavingLAr")->GetValue(j+1);
            ZPosLArEndTrue = truth->GetLeaf("RecoTrackPrimaryParticleTruePositionLeavingLAr")->GetValue(j+2);

            XMomTMSStartTrue = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumEnteringTMS")->GetValue(j);
            YMomTMSStartTrue = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumEnteringTMS")->GetValue(j+1);
            ZMomTMSStartTrue = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumEnteringTMS")->GetValue(j+2);

            XMomLArEndTrue = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumLeavingLAr")->GetValue(j);
            YMomLArEndTrue = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumLeavingLAr")->GetValue(j+1);
            ZMomLArEndTrue = truth->GetLeaf("RecoTrackPrimaryParticleTrueMomentumLeavingLAr")->GetValue(j+2);

            //pulling reco variables
            XPosTMSStartReco = reco->GetLeaf("StartPos")->GetValue(k);
            YPosTMSStartReco = reco->GetLeaf("StartPos")->GetValue(k+1);
            ZPosTMSStartReco = reco->GetLeaf("StartPos")->GetValue(k+2);

            //I think these are intended to be directional cosines (dx, dy, dz) but there is some odd renormalization going on
            XDirTMSStartReco = reco->GetLeaf("StartDirection")->GetValue(k); 
            YDirTMSStartReco = reco->GetLeaf("StartDirection")->GetValue(k+1); 
            ZDirTMSStartReco = reco->GetLeaf("StartDirection")->GetValue(k+2); 

            //for containment
            XEndTrue = truth->GetLeaf("DeathPosition")->GetValue(j);
            YEndTrue = truth->GetLeaf("DeathPosition")->GetValue(j+1);
            ZEndTrue = truth->GetLeaf("DeathPosition")->GetValue(j+2);

            XStartTrue = truth->GetLeaf("BirthPosition")->GetValue(j);
            YStartTrue = truth->GetLeaf("BirthPosition")->GetValue(j+1);
            ZStartTrue = truth->GetLeaf("BirthPosition")->GetValue(j+2);

            float EnergyLeavingLAr, EnergyEnteringTMS;

            EnergyEnteringTMS = truth->GetLeaf("MomentumTMSStart")->GetValue(j+3);
            EnergyLeavingLAr = truth->GetLeaf("MomentumLArEnd")->GetValue(j+3); 

            if(EnergyEnteringTMS < -99999 || EnergyEnteringTMS > 99999){
                continue;
            }

            if(EnergyLeavingLAr < -99999 || EnergyLeavingLAr > 99999){
                continue;
            }

            float deltaE;

            deltaE = EnergyLeavingLAr - EnergyEnteringTMS;
            energyhisto->Fill(deltaE);

            //ensuring normalization

            float TMSRecoDirLen = sqrt((XDirTMSStartReco * XDirTMSStartReco) + (YDirTMSStartReco * YDirTMSStartReco) + (ZDirTMSStartReco * ZDirTMSStartReco));

            /* //for hough lines drawn backwards
            XDirTMSStartReco = (-1) * XDirTMSStartReco / TMSRecoDirLen;
            YDirTMSStartReco = (-1) * YDirTMSStartReco / TMSRecoDirLen;
            ZDirTMSStartReco = (-1) * ZDirTMSStartReco / TMSRecoDirLen;
            */
            
            
            XDirTMSStartReco = XDirTMSStartReco / TMSRecoDirLen;
            YDirTMSStartReco = YDirTMSStartReco / TMSRecoDirLen;
            ZDirTMSStartReco = ZDirTMSStartReco / TMSRecoDirLen;
            

            if(XDirTMSStartReco > 0.99 || XDirTMSStartReco < -0.99){
                //cout << "Y reco failure at event " << reco->GetLeaf("EventNo")->GetValue(0) << " , track number " << (j / 4) + 1 << endl;
                //cout << "X dir: " << XDirTMSStartReco << endl;
                //cout << "Y dir: " << YDirTMSStartReco << endl;
                //cout << "Z dir: " << ZDirTMSStartReco << endl;
                continue;
            }

            //can use to strengthen above cut: && (YDirTMSStartReco < 0.001 && YDirTMSStartReco > -0.001) && (ZDirTMSStartReco < 0.001 && ZDirTMSStartReco > -0.001)

            if(YDirTMSStartReco < 0.001 && YDirTMSStartReco > -0.001){
                cout << "X dir: " << XDirTMSStartReco << endl;
                cout << "Y dir: " << YDirTMSStartReco << endl;
                cout << "Z dir: " << ZDirTMSStartReco << endl;
                continue;
            }

            /*
            //break for events where reco fails
            if(XDirTMSStartReco == 0 && YDirTMSStartReco == 0 && ZDirTMSStartReco == 1){
                cout << "Reco failure at event " << reco->GetLeaf("EventNo")->GetValue(0) << " , track number " << j << endl;
                continue;
            } 
            */

           float ZLength = ZEndTrue - ZPosTMSStartTrue;

           if(ZLength <= 330){
            continue;
           }

            if(XPosTMSStartTrue <= -9999 || YPosTMSStartTrue <= -9999 || XPosLArEndTrue <= -9999 || YPosLArEndTrue <= -9999){ //cuts nonsense
                continue;
            }

            if(YPosTMSStartTrue > 500){ //cuts nonsense
                continue;
            }

            /*if(nHits <= 20){ //cuts tracks with less than 20 hits
                continue;
            }*/

            //containment coords from xiaoyan
           if(XStartTrue > 3700 || XStartTrue < -4500 || YStartTrue > 1000 || YStartTrue < -3200 || ZStartTrue < 4100 || ZStartTrue > 9200){
                continue;
           }

           if(XEndTrue < -3300 || XEndTrue > 3300 || YEndTrue < -2850 || YEndTrue > 160 || ZEndTrue < 11362){ 
            continue;
           }

           //z end containment: 13500 for only thin, 18314 else
           if(ZEndTrue > 18314){
            continue;
           }

            //calculating true starting directions in TMS
            float MomentumLengthTMS = sqrt((XMomTMSStartTrue * XMomTMSStartTrue) + (YMomTMSStartTrue * YMomTMSStartTrue) + (ZMomTMSStartTrue * ZMomTMSStartTrue));
            float dx_tms_true = (XMomTMSStartTrue / MomentumLengthTMS);
			float dy_tms_true = (YMomTMSStartTrue / MomentumLengthTMS);
            float dz_tms_true = (ZMomTMSStartTrue / MomentumLengthTMS);

            //calculating true ending directions in LAr
            float MomentumLengthLAr = sqrt((XMomLArEndTrue * XMomLArEndTrue) + (YMomLArEndTrue * YMomLArEndTrue) + (ZMomLArEndTrue * ZMomLArEndTrue));
            float dx_lar_true = (XMomLArEndTrue / MomentumLengthLAr);
			float dy_lar_true = (YMomLArEndTrue / MomentumLengthLAr);
            float dz_lar_true = (ZMomLArEndTrue / MomentumLengthLAr);

            float xz_angle_true_lar, yz_angle_true_lar, xz_angle_true_tms, yz_angle_true_tms, xz_angle_reco_tms, yz_angle_reco_tms;

            //gives angle between the projection of the relevant vector in the xz or yz plane and the positive z-axis
            xz_angle_true_lar = delta_theta(dx_lar_true, 0, dz_lar_true, 0, 0, 1);
            yz_angle_true_lar = delta_theta(0, dy_lar_true, dz_lar_true, 0, 0, 1);

            if(dx_lar_true < 0){
                xz_angle_true_lar = -1 * xz_angle_true_lar;
            }
            if(dy_lar_true < 0){
                yz_angle_true_lar = -1 * yz_angle_true_lar;
            }

            xz_angle_true_tms = delta_theta(dx_tms_true, 0, dz_tms_true, 0, 0, 1);
            yz_angle_true_tms = delta_theta(0, dy_tms_true, dz_tms_true, 0, 0, 1);

            if(dx_tms_true < 0){
                xz_angle_true_tms = -1 * xz_angle_true_tms;
            }
            if(dy_tms_true < 0){
                yz_angle_true_tms = -1 * yz_angle_true_tms;
            }

            xz_angle_reco_tms = delta_theta(XDirTMSStartReco, 0, ZDirTMSStartReco, 0, 0, 1);
            yz_angle_reco_tms = delta_theta(0, YDirTMSStartReco, ZDirTMSStartReco, 0, 0, 1);

            if(XDirTMSStartReco < 0){
                xz_angle_reco_tms = -1 * xz_angle_reco_tms;
            }
            if(YDirTMSStartReco < 0){
                yz_angle_reco_tms = -1 * yz_angle_reco_tms;
            }

            //filling input histograms
            /*
            x_pos_tms_truth_h->Fill(XPosTMSStartTrue);
            y_pos_tms_truth_h->Fill(YPosTMSStartTrue);
            z_pos_tms_truth_h->Fill(ZPosTMSStartTrue);

            x_pos_lar_truth_h->Fill(XPosLArEndTrue);
            y_pos_lar_truth_h->Fill(YPosLArEndTrue);
            z_pos_lar_truth_h->Fill(ZPosLArEndTrue);

            x_pos_tms_reco_h->Fill(XPosTMSStartReco);
            y_pos_tms_reco_h->Fill(YPosTMSStartReco);
            z_pos_tms_reco_h->Fill(ZPosTMSStartReco);

            dx_tms_true_h->Fill(dx_tms_true);
            dy_tms_true_h->Fill(dy_tms_true);
            dz_tms_true_h->Fill(dz_tms_true);

            dx_lar_true_h->Fill(dx_lar_true);
            dy_lar_true_h->Fill(dy_lar_true);
            dz_lar_true_h->Fill(dz_lar_true);

            x_dir_tms_reco_h->Fill(XDirTMSStartReco);
            y_dir_tms_reco_h->Fill(YDirTMSStartReco);
            z_dir_tms_reco_h->Fill(ZDirTMSStartReco);
            */

            /*
            delta_xz_tms_truereco->Fill((xz_angle_true_tms - xz_angle_reco_tms));
            delta_yz_tms_truereco->Fill((yz_angle_true_tms - yz_angle_reco_tms));

            delta_xz_tmslar_true->Fill((xz_angle_true_tms - xz_angle_true_lar));
            delta_yz_tmslar_true->Fill((yz_angle_true_tms - yz_angle_true_lar));

            delta_xz_tmslar_recotrue->Fill((xz_angle_reco_tms - xz_angle_true_lar));
            delta_yz_tmslar_recotrue->Fill((yz_angle_reco_tms - yz_angle_true_lar));

            //true TMS -> true LAr
            xhisto_tmstolar_true->Fill(delta_x(XPosTMSStartTrue, XPosLArEndTrue, dx_tms_true, dz_tms_true, tmstolar));
            yhisto_tmstolar_true->Fill(delta_y(YPosTMSStartTrue, YPosLArEndTrue, dy_tms_true, dz_tms_true, tmstolar));
            rhisto_tmstolar_true->Fill(delta_r(XPosTMSStartTrue, YPosTMSStartTrue, XPosLArEndTrue, YPosLArEndTrue, dx_tms_true, dy_tms_true, dz_tms_true, tmstolar)); 
            
            //true LAr -> true TMS
            xhisto_lartotms_true->Fill(delta_x(XPosLArEndTrue, XPosTMSStartTrue, dx_lar_true, dz_lar_true, lartotms));
            yhisto_lartotms_true->Fill(delta_y(YPosLArEndTrue, YPosTMSStartTrue, dy_lar_true, dz_lar_true, lartotms));
            rhisto_lartotms_true->Fill(delta_r(XPosLArEndTrue, YPosLArEndTrue, XPosTMSStartTrue, YPosTMSStartTrue, dx_lar_true, dy_lar_true, dz_lar_true, lartotms));
            
            //reco TMS -> true LAr
            xhisto_tmstolar_recototrue->Fill(delta_x(XPosTMSStartReco, XPosLArEndTrue, XDirTMSStartReco, ZDirTMSStartReco, tmstolar));
            yhisto_tmstolar_recototrue->Fill(delta_y(YPosTMSStartReco, YPosLArEndTrue, YDirTMSStartReco, ZDirTMSStartReco, tmstolar));
            rhisto_tmstolar_recototrue->Fill(delta_r(XPosTMSStartReco, YPosTMSStartReco, XPosLArEndTrue, YPosLArEndTrue, XDirTMSStartReco, YDirTMSStartReco, ZDirTMSStartReco, tmstolar));

            //true LAr -> reco TMS
            xhisto_lartotms_truetoreco->Fill(delta_x(XPosLArEndTrue, XPosTMSStartReco, dx_lar_true, dz_lar_true, lartotms));
            yhisto_lartotms_truetoreco->Fill(delta_y(YPosLArEndTrue, YPosTMSStartReco, dy_lar_true, dz_lar_true, lartotms));
            rhisto_lartotms_truetoreco->Fill(delta_r(XPosLArEndTrue, YPosLArEndTrue, XPosTMSStartReco, YPosTMSStartReco, dx_lar_true, dy_lar_true, dz_lar_true, lartotms));
            */
            //delta theta histogram
            //thetahisto_tms_truevreco->Fill(delta_theta(dx_tms_true, dy_tms_true, dz_tms_true, XDirTMSStartReco, YDirTMSStartReco, ZDirTMSStartReco));
            /*
            thetahisto_tmslar_truth->Fill(delta_theta(dx_tms_true, dy_tms_true, dz_tms_true, dx_lar_true, dy_lar_true, dz_lar_true));
            thetahisto_tmslar_recotrue->Fill(delta_theta(XDirTMSStartReco, YDirTMSStartReco, ZDirTMSStartReco, dx_lar_true, dy_lar_true, dz_lar_true));

            deltay_vs_tracklength->Fill(delta_y(YPosTMSStartReco, YPosLArEndTrue, YDirTMSStartReco, ZDirTMSStartReco, tmstolar), ZLength);
            y_dir_tms_reco_vs_tracklength->Fill(YDirTMSStartReco, ZLength);
            delta_theta_tmslar_recotrue_vs_tracklength->Fill(delta_theta(XDirTMSStartReco, YDirTMSStartReco, ZDirTMSStartReco, dx_lar_true, dy_lar_true, dz_lar_true), ZLength);

            theta_xz_tmsreco_vs_theta_xz_lartruth->Fill(xz_angle_true_lar, xz_angle_reco_tms);
            theta_yz_tmsreco_vs_theta_yz_lartruth->Fill(yz_angle_true_lar, yz_angle_reco_tms);
            */
        }
    }

    //draw & format input histograms
    /*
    x_pos_tms_truth_h->SetTitle("True Starting Position in X in TMS");
    x_pos_tms_truth_h->Draw();
    std::string printname = "plots/" + filename + "/1d/tms_xpos_truth.png";
    c1->Print(printname.c_str());

    y_pos_tms_truth_h->SetTitle("True Starting Position in Y in TMS");
    y_pos_tms_truth_h->Draw();
    printname = "plots/" + filename + "/1d/tms_ypos_truth.png";
    c1->Print(printname.c_str());

    z_pos_tms_truth_h->SetTitle("True Starting Position in Z in TMS");
    z_pos_tms_truth_h->Draw();
    printname = "plots/" + filename + "/1d/tms_zpos_truth.png";
    c1->Print(printname.c_str());

    x_pos_lar_truth_h->SetTitle("True Ending Position in X in LAr");
    x_pos_lar_truth_h->Draw();
    printname = "plots/" + filename + "/1d/lar_xpos_truth.png";
    c1->Print(printname.c_str());

    y_pos_lar_truth_h->SetTitle("True Ending Position in Y in LAr");
    y_pos_lar_truth_h->Draw();
    printname = "plots/" + filename + "/1d/lar_ypos_truth.png";
    c1->Print(printname.c_str());

    z_pos_lar_truth_h->SetTitle("True Ending Position in Z in LAr");
    z_pos_lar_truth_h->Draw();
    printname = "plots/" + filename + "/1d/lar_zpos_truth.png";
    c1->Print(printname.c_str());

    x_pos_tms_reco_h->SetTitle("Reconstructed Starting Position in X in TMS");
    x_pos_tms_reco_h->Draw();
    printname = "plots/" + filename + "/1d/tmsreco_xpos.png";
    c1->Print(printname.c_str());

    y_pos_tms_reco_h->SetTitle("Reconstructed Starting Position in Y in TMS");
    y_pos_tms_reco_h->Draw();
    printname = "plots/" + filename + "/1d/tmsreco_ypos.png";
    c1->Print(printname.c_str());

    z_pos_tms_reco_h->SetTitle("Reconstructed Starting Position in Z in TMS");
    z_pos_tms_reco_h->Draw();
    printname = "plots/" + filename + "/1d/tmsreco_zpos.png";
    c1->Print(printname.c_str());

    dx_tms_true_h->SetTitle("X Momentum at TMS Start / |P_TMS|");
    dx_tms_true_h->Draw();
    printname = "plots/" + filename + "/1d/tms_xdir_true.png";
    c1->Print(printname.c_str());

    dy_tms_true_h->SetTitle("Y Momentum at TMS Start / |P_TMS|");
    dy_tms_true_h->Draw();
    printname = "plots/" + filename + "/1d/tms_ydir_true.png";
    c1->Print(printname.c_str());

    dz_tms_true_h->SetTitle("Z Momentum at TMS Start / |P_TMS|");
    dz_tms_true_h->Draw();
    printname = "plots/" + filename + "/1d/tms_zdir_true.png";
    c1->Print(printname.c_str());

    dx_lar_true_h->SetTitle("X Momentum at LAr End / |P_LAr|");
    dx_lar_true_h->Draw();
    printname = "plots/" + filename + "/1d/lar_xdir_true.png";
    c1->Print(printname.c_str());

    dy_lar_true_h->SetTitle("Y Momentum at LAr End / |P_LAr|");
    dy_lar_true_h->Draw();
    printname = "plots/" + filename + "/1d/lar_ydir_true.png";
    c1->Print(printname.c_str());

    dz_lar_true_h->SetTitle("Z Momentum at LAr End / |P_LAr|");
    dz_lar_true_h->Draw();
    printname = "plots/" + filename + "/1d/lar_zdir_true.png";
    c1->Print(printname.c_str());

    x_dir_tms_reco_h->SetTitle("StartDirection (X), TMS Reco");
    x_dir_tms_reco_h->Draw();
    printname = "plots/" + filename + "/1d/tms_xdir_reco.png";
    c1->Print(printname.c_str());

    y_dir_tms_reco_h->SetTitle("StartDirection (Y), TMS Reco");
    y_dir_tms_reco_h->Draw();
    printname = "plots/" + filename + "/1d/tms_ydir_reco.png";
    c1->Print(printname.c_str());

    z_dir_tms_reco_h->SetTitle("StartDirection (Z), TMS Reco");
    z_dir_tms_reco_h->Draw();
    printname = "plots/" + filename + "/1d/tms_zdir_reco.png";
    c1->Print(printname.c_str());
    */
    /*
    delta_xz_tms_truereco->SetTitle("Angle between True TMS Dir. and Reco TMS Dir. in the XZ Plane");
    delta_xz_tms_truereco->GetXaxis()->SetTitle("Degrees");
    delta_xz_tms_truereco->Draw();
    printname = "plots/" + filename + "/delta_xz_tms_truereco.png";
    c1->Print(printname.c_str());

    delta_yz_tms_truereco->SetTitle("Angle between True TMS Dir. and Reco TMS Dir. in the YZ Plane");
    delta_yz_tms_truereco->GetXaxis()->SetTitle("Degrees");
    delta_yz_tms_truereco->Draw();
    printname = "plots/" + filename + "/delta_yz_tms_truereco.png";
    c1->Print(printname.c_str());

    delta_xz_tmslar_true->SetTitle("Angle between True TMS Dir. and True LAr Dir. in the XZ Plane");
    delta_xz_tmslar_true->GetXaxis()->SetTitle("Degrees");
    delta_xz_tmslar_true->Draw();
    printname = "plots/" + filename + "/delta_xz_tmslar_true.png";
    c1->Print(printname.c_str());

    delta_yz_tmslar_true->SetTitle("Angle between True TMS Dir. and True LAr Dir. in the YZ Plane");
    delta_yz_tmslar_true->GetXaxis()->SetTitle("Degrees");
    delta_yz_tmslar_true->Draw();
    printname = "plots/" + filename + "/delta_yz_tmslar_true.png";
    c1->Print(printname.c_str());

    delta_xz_tmslar_recotrue->SetTitle("Angle between Reco TMS Dir. and True LAr Dir. in the XZ Plane");
    delta_xz_tmslar_recotrue->GetXaxis()->SetTitle("Degrees");
    delta_xz_tmslar_recotrue->Draw();
    printname = "plots/" + filename + "/delta_xz_tmslar_recotrue.png";
    c1->Print(printname.c_str());

    delta_yz_tmslar_recotrue->SetTitle("Angle between Reco TMS Dir. and True LAr Dir. in the YZ Plane");
    delta_yz_tmslar_recotrue->GetXaxis()->SetTitle("Degrees");
    delta_yz_tmslar_recotrue->Draw();
    printname = "plots/" + filename + "/delta_yz_tmslar_recotrue.png";
    c1->Print(printname.c_str());

    //draw, fit, & format delta histograms
    xhisto_tmstolar_true->SetTitle("\\Delta X");
    xhisto_tmstolar_true->GetXaxis()->SetTitleOffset(1.25);
    xhisto_tmstolar_true->GetXaxis()->SetTitle("Extrapolated True TMS Start Pos - True LAr End Pos (mm)");
    xhisto_tmstolar_true->GetYaxis()->SetTitleOffset(1.25);
    xhisto_tmstolar_true->GetYaxis()->SetTitle("N Tracks");
    xhisto_tmstolar_true->Fit("gaus");
    TF1 *fitobj = (TF1*)xhisto_tmstolar_true->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    fitinfo->SetTextSizePixels(35);
    TText* head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    xhisto_tmstolar_true->GetYaxis()->SetRangeUser(xhisto_tmstolar_true->GetYaxis()->GetXmin(), xhisto_tmstolar_true->GetMaximum()*1.25);
    dunestyle::CenterTitles(xhisto_tmstolar_true);
    xhisto_tmstolar_true->Draw();
    fitinfo->Draw(); 
    printname = "plots/" + filename + "/deltax_tmstolar_truth.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    yhisto_tmstolar_true->SetTitle("\\Delta Y");
    yhisto_tmstolar_true->GetXaxis()->SetTitleOffset(1.25);
    yhisto_tmstolar_true->GetXaxis()->SetTitle("Extrapolated True TMS Start Pos - True LAr End Pos (mm)");
    yhisto_tmstolar_true->GetYaxis()->SetTitleOffset(1.25);
    yhisto_tmstolar_true->GetYaxis()->SetTitle("N Tracks");
    yhisto_tmstolar_true->Fit("gaus");
    fitobj = (TF1*)yhisto_tmstolar_true->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    fitinfo->SetTextSizePixels(35);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    yhisto_tmstolar_true->GetYaxis()->SetRangeUser(yhisto_tmstolar_true->GetYaxis()->GetXmin(), yhisto_tmstolar_true->GetMaximum()*1.25);
    dunestyle::CenterTitles(yhisto_tmstolar_true);
    yhisto_tmstolar_true->Draw();
    fitinfo->Draw();
    printname = "plots/" + filename + "/deltay_tmstolar_truth.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    rhisto_tmstolar_true->SetTitle("\\Delta R");
    rhisto_tmstolar_true->GetXaxis()->SetTitle("Extrapolated True TMS Start Pos - True LAr End Pos (mm)");
    rhisto_tmstolar_true->Draw();
    printname = "plots/" + filename + "/deltar_tmstolar_truth.png";
    c1->Print(printname.c_str());

    xhisto_lartotms_true->SetTitle("\\Delta X");
    xhisto_lartotms_true->GetXaxis()->SetTitle("Extrapolated True LAr End Pos - True TMS Start Pos (mm)");
    xhisto_lartotms_true->GetXaxis()->SetTitleOffset(1.25);
    xhisto_lartotms_true->GetYaxis()->SetTitle("N Tracks");
    xhisto_lartotms_true->GetYaxis()->SetTitleOffset(1.25);
    xhisto_lartotms_true->Fit("gaus");
    fitobj = (TF1*)xhisto_lartotms_true->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    fitinfo->SetTextSizePixels(35);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    xhisto_lartotms_true->GetYaxis()->SetRangeUser(xhisto_lartotms_true->GetYaxis()->GetXmin(), xhisto_lartotms_true->GetMaximum()*1.25);
    dunestyle::CenterTitles(xhisto_lartotms_true);
    xhisto_lartotms_true->Draw();
    fitinfo->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/deltax_lartotms_truth.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    yhisto_lartotms_true->SetTitle("\\Delta Y");
    yhisto_lartotms_true->GetXaxis()->SetTitle("Extrapolated True LAr End Pos - True TMS Start Pos (mm)");
    yhisto_lartotms_true->GetXaxis()->SetTitleOffset(1.25);
    yhisto_lartotms_true->GetYaxis()->SetTitle("N Tracks");
    yhisto_lartotms_true->GetYaxis()->SetTitleOffset(1.25);
    yhisto_lartotms_true->Fit("gaus");
    fitobj = (TF1*)yhisto_lartotms_true->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    fitinfo->SetTextSizePixels(35);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    yhisto_lartotms_true->GetYaxis()->SetRangeUser(yhisto_lartotms_true->GetYaxis()->GetXmin(), yhisto_lartotms_true->GetMaximum()*1.25);
    dunestyle::CenterTitles(yhisto_lartotms_true);
    yhisto_lartotms_true->Draw();
    fitinfo->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/deltay_lartotms_truth.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    rhisto_lartotms_true->SetTitle("\\Delta R");
    rhisto_lartotms_true->GetXaxis()->SetTitle("Extrapolated True LAr End Pos - True TMS Start Pos (mm)");
    rhisto_lartotms_true->Draw();
    printname = "plots/" + filename + "/deltar_lartotms_truth.png";
    c1->Print(printname.c_str());

    xhisto_tmstolar_recototrue->SetTitle("\\Delta X");
    xhisto_tmstolar_recototrue->GetXaxis()->SetTitle("Extrapolated Reco TMS Start Pos - True LAr End Pos (mm)");
    xhisto_tmstolar_recototrue->GetYaxis()->SetTitle("N Tracks");
    xhisto_tmstolar_recototrue->Fit("gaus");
    fitobj = (TF1*)xhisto_tmstolar_recototrue->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    xhisto_tmstolar_recototrue->GetYaxis()->SetRangeUser(xhisto_tmstolar_recototrue->GetYaxis()->GetXmin(), xhisto_tmstolar_recototrue->GetMaximum()*1.25);
    dunestyle::CenterTitles(xhisto_tmstolar_recototrue);
    xhisto_tmstolar_recototrue->Draw();
    fitinfo->Draw();
    printname = "plots/" + filename + "/deltax_tmstolar_recototrue.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    yhisto_tmstolar_recototrue->SetTitle("\\Delta Y");
    yhisto_tmstolar_recototrue->GetXaxis()->SetTitle("Extrapolated Reco TMS Start Pos - True LAr End Pos (mm)");
    yhisto_tmstolar_recototrue->GetYaxis()->SetTitle("N Tracks");
    yhisto_tmstolar_recototrue->Fit("gaus");
    fitobj = (TF1*)yhisto_tmstolar_recototrue->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    yhisto_tmstolar_recototrue->GetYaxis()->SetRangeUser(yhisto_tmstolar_recototrue->GetYaxis()->GetXmin(), yhisto_tmstolar_recototrue->GetMaximum()*1.25);
    dunestyle::CenterTitles(yhisto_tmstolar_recototrue);
    yhisto_tmstolar_recototrue->Draw();
    fitinfo->Draw();
    printname = "plots/" + filename + "/deltay_tmstolar_recototrue.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    rhisto_tmstolar_recototrue->SetTitle("\\Delta R");
    rhisto_tmstolar_recototrue->GetXaxis()->SetTitle("Extrapolated Reco TMS Start Pos - True LAr End Pos (mm)");
    rhisto_tmstolar_recototrue->Draw();
    printname = "plots/" + filename + "/deltar_tmstolar_recototrue.png";
    c1->Print(printname.c_str());

    xhisto_lartotms_truetoreco->SetTitle("\\Delta X");
    xhisto_lartotms_truetoreco->GetXaxis()->SetTitle("Extrapolated True LAr End Pos - Reco TMS Start Pos (mm)");
    xhisto_lartotms_truetoreco->GetYaxis()->SetTitleOffset(1.25);
    xhisto_lartotms_truetoreco->GetYaxis()->SetTitle("N Tracks");
    xhisto_lartotms_truetoreco->GetYaxis()->SetTitleOffset(1.25);
    xhisto_lartotms_truetoreco->Fit("gaus");
    fitobj = (TF1*)xhisto_lartotms_truetoreco->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    fitinfo->SetTextSizePixels(35);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    xhisto_lartotms_truetoreco->GetYaxis()->SetRangeUser(xhisto_lartotms_truetoreco->GetYaxis()->GetXmin(), xhisto_lartotms_truetoreco->GetMaximum()*1.25);
    dunestyle::CenterTitles(xhisto_lartotms_truetoreco);
    xhisto_lartotms_truetoreco->Draw();
    fitinfo->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/deltax_lartotms_truetoreco.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    yhisto_lartotms_truetoreco->SetTitle("\\Delta Y");
    yhisto_lartotms_truetoreco->GetXaxis()->SetTitle("Extrapolated True LAr End Pos - Reco TMS Start Pos (mm)");
    yhisto_lartotms_truetoreco->GetYaxis()->SetTitleOffset(1.25);
    yhisto_lartotms_truetoreco->GetYaxis()->SetTitle("N Tracks");
    yhisto_lartotms_truetoreco->GetYaxis()->SetTitleOffset(1.25);
    yhisto_lartotms_truetoreco->Fit("gaus");
    fitobj = (TF1*)yhisto_lartotms_truetoreco->GetListOfFunctions()->FindObject("gaus");
    mean = fitobj->GetParameter(1);
    sigma = fitobj->GetParameter(2);
    fitinfo->SetTextSizePixels(35);
    head = fitinfo->AddText("Fit Parameters:");
    head->SetTextFont(62);
    fitinfo->AddText(Form("Mean: %g", mean));
    fitinfo->AddText(Form("\\sigma: %g", sigma));
    yhisto_lartotms_truetoreco->GetYaxis()->SetRangeUser(yhisto_lartotms_truetoreco->GetYaxis()->GetXmin(), yhisto_lartotms_truetoreco->GetMaximum()*1.25);
    dunestyle::CenterTitles(yhisto_lartotms_truetoreco);
    yhisto_lartotms_truetoreco->Draw();
    fitinfo->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/deltay_lartotms_truetoreco.png";
    c1->Print(printname.c_str());
    fitinfo->Clear();

    rhisto_lartotms_truetoreco->SetTitle("\\Delta R");
    rhisto_lartotms_truetoreco->GetXaxis()->SetTitle("Extrapolated True LAr End Pos - Reco TMS Start Pos (mm)");
    rhisto_lartotms_truetoreco->Draw();
    printname = "plots/" + filename + "/deltar_lartotms_truetoreco.png";
    c1->Print(printname.c_str());

    //draw, fit, & print delta theta histograms
    
    thetahisto_tms_truevreco->SetTitle("\\Delta\\theta");
    thetahisto_tms_truevreco->GetXaxis()->SetTitle("Angle Between True and Reco TMS Start Directions (degrees)");
    thetahisto_tms_truevreco->GetYaxis()->SetTitle("N Tracks");
    thetahisto_tms_truevreco->GetXaxis()->SetTitleSize(0.035);
    thetahisto_tms_truevreco->GetXaxis()->SetTitleOffset(1.2);
    thetahisto_tms_truevreco->GetYaxis()->SetTitleSize(0.035);
    thetahisto_tms_truevreco->GetYaxis()->SetTitleOffset(1.7);
    dunestyle::CenterTitles(thetahisto_tms_truevreco);
    thetahisto_tms_truevreco->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/thetahisto_tms_truevreco.png";
    c1->Print(printname.c_str());
    
    
    thetahisto_tmslar_truth->SetTitle("\\Delta\\theta");
    thetahisto_tmslar_truth->GetYaxis()->SetTitle("N Tracks");
    thetahisto_tmslar_truth->GetXaxis()->SetTitle("Angle Between True TMS Start Direction/True LAr End Direction (degrees)");
    thetahisto_tmslar_truth->GetXaxis()->SetTitleSize(0.035);
    thetahisto_tmslar_truth->GetXaxis()->SetTitleOffset(1.2);
    thetahisto_tmslar_truth->GetYaxis()->SetTitleSize(0.035);
    thetahisto_tmslar_truth->GetYaxis()->SetTitleOffset(1.7);
    dunestyle::CenterTitles(thetahisto_tmslar_truth);
    thetahisto_tmslar_truth->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/thetahisto_tmslar_truth.png";
    c1->Print(printname.c_str());

    thetahisto_tmslar_recotrue->SetTitle("\\Delta\\theta");
    thetahisto_tmslar_recotrue->GetXaxis()->SetTitle("Angle Between Reco TMS Start Direction/True LAr End Direction (degrees)");
    thetahisto_tmslar_recotrue->GetYaxis()->SetTitle("N Tracks");
    thetahisto_tmslar_recotrue->GetXaxis()->SetTitleSize(0.035);
    thetahisto_tmslar_recotrue->GetXaxis()->SetTitleOffset(1.2);
    thetahisto_tmslar_recotrue->GetYaxis()->SetTitleSize(0.035);
    thetahisto_tmslar_recotrue->GetYaxis()->SetTitleOffset(1.7);
    dunestyle::CenterTitles(thetahisto_tmslar_recotrue);
    thetahisto_tmslar_recotrue->Draw();
    dunestyle::Simulation();
    printname = "plots/" + filename + "/thetahisto_tmslar_recotrue.png";
    c1->Print(printname.c_str());

    deltay_vs_tracklength->SetTitle("Z Length vs. \\Delta Y");
    deltay_vs_tracklength->GetXaxis()->SetTitle("Extrapolated True LAr End Position - Reco TMS Start Position (mm)");
    deltay_vs_tracklength->GetYaxis()->SetTitle("Length in TMS in Z (mm)");
    deltay_vs_tracklength->GetXaxis()->SetTitleSize(0.035);
    deltay_vs_tracklength->GetXaxis()->SetTitleOffset(1.2);
    deltay_vs_tracklength->GetYaxis()->SetTitleSize(0.035);
    deltay_vs_tracklength->GetYaxis()->SetTitleOffset(1.7);
    dunestyle::CenterTitles(deltay_vs_tracklength);
    deltay_vs_tracklength->Draw("COLZ");
    dunestyle::Simulation();
    printname = "plots/" + filename + "/deltay_vs_tracklength.png";
    c1->Print(printname.c_str());

    y_dir_tms_reco_vs_tracklength->SetTitle("Z Length vs. Reco TMS Start Y Dir.");
    y_dir_tms_reco_vs_tracklength->GetXaxis()->SetTitle("Reco TMS Start Direction (Y)");
    y_dir_tms_reco_vs_tracklength->GetYaxis()->SetTitle("Length in TMS in Z (mm)");
    y_dir_tms_reco_vs_tracklength->Draw("COLZ");
    dunestyle::Simulation();
    printname = "plots/" + filename + "/y_dir_tms_reco_vs_tracklength.png";
    c1->Print(printname.c_str());

    TLine *l = new TLine(-30.0, -30.0, 30.0, 30.0);
    l->SetLineColor(2);

    theta_xz_tmsreco_vs_theta_xz_lartruth->SetTitle("TMS Reco \\theta_{xz} vs. ND-LAr Truth \\theta_{xz}");
    theta_xz_tmsreco_vs_theta_xz_lartruth->GetXaxis()->SetTitle("ND-LAr Truth \\theta_{xz} (degrees)");
    theta_xz_tmsreco_vs_theta_xz_lartruth->GetYaxis()->SetTitle("TMS Reco \\theta_{xz} (degrees)");
    dunestyle::CenterTitles(theta_xz_tmsreco_vs_theta_xz_lartruth);
    theta_xz_tmsreco_vs_theta_xz_lartruth->Draw("COLZ");
    l->Draw("SAME");
    dunestyle::Simulation();
    printname = "plots/" + filename + "/theta_xz_tmsreco_vs_theta_xz_lartruth.png";
    c1->Print(printname.c_str());

    theta_yz_tmsreco_vs_theta_yz_lartruth->SetTitle("TMS Reco \\theta_{yz} vs. ND-LAr Truth \\theta_{yz}");
    theta_yz_tmsreco_vs_theta_yz_lartruth->GetXaxis()->SetTitle("ND-LAr Truth \\theta_{yz} (degrees)");
    theta_yz_tmsreco_vs_theta_yz_lartruth->GetYaxis()->SetTitle("TMS Reco \\theta_{yz} (degrees)");
    dunestyle::CenterTitles(theta_yz_tmsreco_vs_theta_yz_lartruth);
    theta_yz_tmsreco_vs_theta_yz_lartruth->Draw("COLZ");
    l->Draw("SAME");
    dunestyle::Simulation();
    printname = "plots/" + filename + "/theta_yz_tmsreco_vs_theta_yz_lartruth.png";
    c1->Print(printname.c_str());

    delta_theta_tmslar_recotrue_vs_tracklength->SetTitle("Z Length in TMS vs. \\Delta\\theta");
    delta_theta_tmslar_recotrue_vs_tracklength->GetYaxis()->SetTitle("Length in TMS in Z (mm)");
    delta_theta_tmslar_recotrue_vs_tracklength->GetXaxis()->SetTitle("Angle Between Reco TMS Start Direction/True LAr End Direction (degrees)");
    delta_theta_tmslar_recotrue_vs_tracklength->GetXaxis()->SetTitleSize(0.035);
    delta_theta_tmslar_recotrue_vs_tracklength->GetXaxis()->SetTitleOffset(1.2);
    delta_theta_tmslar_recotrue_vs_tracklength->GetYaxis()->SetTitleSize(0.035);
    delta_theta_tmslar_recotrue_vs_tracklength->GetYaxis()->SetTitleOffset(1.9);
    dunestyle::CenterTitles(delta_theta_tmslar_recotrue_vs_tracklength);
    delta_theta_tmslar_recotrue_vs_tracklength->Draw("COLZ");
    dunestyle::Simulation();
    printname = "plots/" + filename + "/delta_theta_tmslar_recotrue_vs_tracklength.png";
    c1->Print(printname.c_str());
    */

    energyhisto->SetTitle("Muon Energy Leaving LAr - Muon Energy Entering TMS");
    energyhisto->Draw();
    std::string printname = "plots/" + filename + "energyhisto.png";
    c1->Print(printname.c_str());
}