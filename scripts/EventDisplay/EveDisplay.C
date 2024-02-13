#include "TGeoManager.h"
#include "TEveManager.h"
#include "TEveTrackPropagator.h"


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
    //gEve->AddGlobalElement(TMSTopNode);
    gEve->GetGlobalScene()->AddElement(TMSTopNode);

    //TEveTrack* track = new TEveTrack();
    //TEveStraightLineSet* tracks = new TEveStraightLineSet("loins");
    //std::cout << tracks << std::endl;
    //tracks->Print();
    //tracks->SetLine(0,0,0,0,0,1,1);


    TEveElement* evt = gEve->GetEventScene();

    auto prop = new TEveTrackPropagator();
    prop->SetMaxR(400);
    prop->SetMaxZ(400);
    prop->SetMagField(0,0,0);
    prop->SetMaxOrbs(3);
    auto tHolder = new TEveCompound("traeks");
    TParticle p;
    p.SetProductionVertex(0,0,-300,1);


    float startPos[3];
    float endPos[3];
    float direction[3];
    float length;

    double Theta;
    double Phi;
    double momScale = 1.0;

    TFile* RecoFile = new TFile("out.root");
    TTree* Reco = RecoFile->Get("Reco_Tree");

    Reco->Branch("StartPos", startPos, "StartPos/F");
    Reco->Branch("EndPos", endPos, "EndPos/F");
    Reco->Branch("Direction", direction, "Direction/F");
    Reco->Branch("Length", &length, "Length/F");

//    for (int i=0; i<3; i++)
//        trackLen[i] = endPos[i] - startPos[i];


    p.SetPolarPhi(Phi);
    p.SetPolarTheta(Theta);
    p.SetMomentum(momScale*direction[0],momScale*direction[1],momScale*direction[2],1);

    p.SetLabel("mu");
    p.SetPdgCode(12);
    //p.SetCharge(1);
    auto t = new TEveTrack(&p, 1, prop);
    t->SetName("AXIONS CONFIRMED");
    t->SetMainColor(kOrange);
    t->SetLineWidth(4);
    t->MakeTrack();
    tHolder->AddElement(t);
    std::cout << prop->GetTrackLength() << std::endl;


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

    TEveViewer *ev = gEve->GetDefaultViewer();
    TGLViewer* gv = ev->GetGLViewer();
    gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, nullptr);

    evt->AddElement(tHolder);

    auto vnt = new TEveEventManager("vent", "fuck TEve and fuck ROOT in general");
    gEve->AddEvent(vnt);
    vnt->NextEvent();

    gEve->Redraw3D(kTRUE);
    tHolder->Paint();


    return 0;
}




