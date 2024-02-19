#include "TGeoManager.h"
#include "TEveManager.h"
#include "TEveTrackPropagator.h"


void DrawSpill(int iEvt = 0)
{
    0x90;
    return;
}

void DropSpill()
{
    gEve->GetViewers()->DeleteAnnotations();   
    gEve->GetCurrentEvent()->DestroyElements();
}

bool GoToSpill()
{
  DropSpill();
  gEve->Redraw3D(kFALSE, kTRUE);
  return true;
}

int EveDisplay()
{
    TEveManager::Create();

    gEve->RegisterGeometryAlias("World", "TMS.root");

    TGeoManager* Geom  = gEve->GetGeometryByAlias("World");//->Import("TMS.gdml");

    TGeoNode* WorldNode = gGeoManager->GetTopVolume()->FindNode("rockBox_lv_0"); // 1 dottir
    TGeoNode* HallNode  = WorldNode->GetDaughter(0);// 8 dottirs
    TGeoNode* TMSNode = HallNode->GetDaughter(1); // TMS is index 1 (2nd only to LAr, SAND in shambles)

    TGeoNode* SteelNodeThinS = TMSNode->GetDaughter(0)->GetDaughter(0); // side panels
    TGeoNode* SteelNodeThinM = TMSNode->GetDaughter(0)->GetDaughter(2); // mid panels
    TGeoNode* SteelNodeThicS = TMSNode->GetDaughter(40)->GetDaughter(0); // steel thicc
    TGeoNode* SteelNodeThicM = TMSNode->GetDaughter(40)->GetDaughter(2);

    TGeoNode* ScintNodeR = TMSNode->GetDaughter(99);//->GetDaughter(0);
    TGeoNode* ScintNodeL = TMSNode->GetDaughter(110);//->GetDaughter(0); // steel thicc

    SteelNodeThinS->GetVolume()->SetOption("box");
    SteelNodeThicM->GetVolume()->SetOption("box");
    SteelNodeThinS->GetVolume()->SetOption("box");
    SteelNodeThicM->GetVolume()->SetOption("box");
    SteelNodeThinS->GetVolume()->SetLineColor(kCyan);
    SteelNodeThinM->GetVolume()->SetLineColor(kCyan);
    SteelNodeThicS->GetVolume()->SetLineColor(kCyan);
    SteelNodeThicM->GetVolume()->SetLineColor(kCyan);
    SteelNodeThinS->GetVolume()->SetTransparency(98);
    SteelNodeThinM->GetVolume()->SetTransparency(98);
    SteelNodeThicS->GetVolume()->SetTransparency(98);
    SteelNodeThicM->GetVolume()->SetTransparency(98);

    //ScintNodeL->GetDaughter(0)->GetVolume()->VisibleDaughters(false);
    ScintNodeR->GetDaughter(0)->GetVolume()->VisibleDaughters(false);

    // same objects for L, R as they're copies of one box; setting one sets both
    ScintNodeR->GetDaughter(0)->GetVolume()->SetOption("w");
    ScintNodeR->GetDaughter(0)->GetVolume()->SetLineColor(kPink);
    ScintNodeR->GetDaughter(0)->GetVolume()->SetTransparency(90);

    TEveGeoTopNode* TMSTopNode= new TEveGeoTopNode(gGeoManager, TMSNode);
    TEveViewer *ev = gEve->GetDefaultViewer();
    TGLViewer* gv = ev->GetGLViewer();
    gEve->GetGlobalScene()->AddElement(TMSTopNode);
    //gEve->AddGlobalElement(TMSTopNode);

    //TEveTrack* track = new TEveTrack();
    //TEveStraightLineSet* tracks = new TEveStraightLineSet("loins");
    //std::cout << tracks << std::endl;
    //tracks->Print();
    //tracks->SetLine(0,0,0,0,0,1,1);



    int nTracksMax = 50; // Hardcoded max tracks

    int nTracks; // tracks in the spill/event
    float startPos[nTracksMax][3];
    float endPos[nTracksMax][3];
    float direction[nTracksMax][3];
    float length;

    double Theta;
    double Phi;
    double momScale = 1.0e3;

    TFile* RecoFile = new TFile("out.root");
    TTree* Reco = (TTree*) RecoFile->Get("Reco_Tree");
    int nEntries = Reco->GetEntries();

    Reco->SetBranchAddress("nTracks", &nTracks);
    Reco->SetBranchAddress("StartPos", startPos);
    Reco->SetBranchAddress("EndPos", endPos);
    Reco->SetBranchAddress("Direction", direction);
    Reco->SetBranchAddress("Length", &length);

    TEveElement* evt = gEve->GetEventScene();

    TEveCompound* tHolder = new TEveCompound("traeks");

    TEveTrack* t;
    TParticle* p;
    //TEveTrackPropagator* prop;
    TEveTrackPropagator* prop = new TEveTrackPropagator();

    for (int i=0; i<nEntries; i++)
    {
        Reco->GetEntry(i);
        if (nTracks <= 0)
        { 
            std::cout << "no tracks in event " << i << ", skipping..." << std::endl;
            continue;
        }

        //TEveTrackPropagator* prop = new TEveTrackPropagator();
        //prop = new TEveTrackPropagator();
        prop->SetMaxR(352);
        prop->SetMaxZ(352);
        prop->SetMagField(0,0,0);
        prop->SetMaxOrbs(3);
        //TParticle* p = new TParticle;
        p = new TParticle;

        std::cout << "startPos: " << startPos[0][0]/10. << " " << startPos[0][1]/10. << " " << startPos[0][2]/10. << std::endl;
        std::cout << "endPos: " << endPos[0][0]/10. << " " << endPos[0][1]/10. << " " << endPos[0][2]/10. << std::endl;
        //std::cout << "dirPos: " << endPos[0][0] - startPos[0][0] << " " << endPos[0][1] - startPos[0][1] << " " << endPos[0][2] - startPos[0][2] << std::endl;
        //std::cout << "dir: " << direction[0][0] << " " << direction[0][1] << " " << direction[0][2] << std::endl;
        //std::cout << "length: " << length << std::endl;
        p->SetProductionVertex(endPos[0][0]/10.,endPos[0][1]/10. + 390.0 ,endPos[0][2]/10. -1485.0 ,1);
        //p->SetProductionVertex(startPos[0][0]/10.,startPos[0][1]/10. +350.0 ,startPos[0][2]/10. -1485.0 ,1);


        p->SetPolarPhi(0);
        p->SetPolarTheta(0);
        p->SetMomentum(-1*momScale*direction[0][0],-1*momScale*direction[0][1],-1*momScale*direction[0][2],1);
        //p->SetMomentum(momScale*direction[0][0],momScale*direction[0][1],momScale*direction[0][2],1);

        //p.SetLabel("mu");
        p->SetPdgCode(12);
        //p.SetCharge(1);
        //TEveTrack* t = new TEveTrack(p, 1, prop);
        t = new TEveTrack(p, 1, prop);
        t->SetName("AXIONS CONFIRMED");
        t->SetMainColor(kOrange);
        t->SetLineWidth(4);
        t->MakeTrack();
        tHolder->AddElement(t);
    }

    /*
    TEveTrackList* tracklist = new TEveTrackList("shit tracks");
    TEveTrackPropagator* tprop = new tracklist->GetPropagator();
    tprop->SetStepper(TEveTrackPropagator::kRungeKutta);
    tracklist->SetMainColor(6);
    tracklist->SetMarkerColor(kYellow);
    tracklist->SetMarkerStyle(4);
    tracklist->SetMarkerSize(0.5);
    TEveTrack* track = new TEveTrack();
    */

    //TEveRecTrack
    //TEvePointSet()


    evt->AddElement(tHolder);
    auto vnt = new TEveEventManager("vent", "fuck TEve and fuck ROOT in general");
    gEve->AddEvent(vnt);
    vnt->NextEvent();

    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, nullptr);
    gEve->Redraw3D(kTRUE);
    tHolder->Paint();
    vnt->NextEvent();


    return 0;
}




