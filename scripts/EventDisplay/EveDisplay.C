#include "TGeoManager.h"
#include "TEveManager.h"

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
    SteelNodeThinS->GetVolume()->SetTransparency(0);
    SteelNodeThinM->GetVolume()->SetTransparency(0);
    SteelNodeThicS->GetVolume()->SetTransparency(0);
    SteelNodeThicM->GetVolume()->SetTransparency(0);

    //ScintNodeL->GetDaughter(0)->GetVolume()->VisibleDaughters(false);
    ScintNodeR->GetDaughter(0)->GetVolume()->VisibleDaughters(false);

    // same objects for L, R as they're copies of one box; setting one sets both
    ScintNodeR->GetDaughter(0)->GetVolume()->SetOption("w");
    ScintNodeR->GetDaughter(0)->GetVolume()->SetLineColor(kPink);
    ScintNodeR->GetDaughter(0)->GetVolume()->SetTransparency(0);

    TEveGeoTopNode* TMSTopNode= new TEveGeoTopNode(gGeoManager, TMSNode);
    gEve->AddGlobalElement(TMSTopNode);

    //TEveTrack* track = new TEveTrack();
    //TEveStraightLineSet* tracks = new TEveStraightLineSet("loins");
    //std::cout << tracks << std::endl;
    //tracks->Print();
    //tracks->SetLine(0,0,0,0,0,1,1);



    gEve->Redraw3D(kTRUE);
    //tracks->Paint();

    return 0;
}
