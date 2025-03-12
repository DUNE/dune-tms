// Copied over from $ROOTSYS/tutorials/geom/geomDO.C
//
// View the world geometry using the TGeoManager
// Works fine in ROOT6, but not ROOT5
// Intended to be played around with to make people comfortable with the geometry rather than treat the code like gospel, hence lots of code is commented out; these are intended as examples
//
// There are many volumes to draw here, so check PrintAllVolumes() output and select the specific volume that you want to draw, if you're so inclined
//
// To run this, do 
// root view_geom.cpp
// or if you want to feed a specific geometry do 
// root 'view_geom.cpp("geometry_file.gdml")'

void RecursiveInvisible(TGeoVolume *vol);
void RecursiveVisible(TGeoVolume *vol);
void RecursiveTransparency(TGeoVolume *vol, Int_t transp);
void PrintRecursive(TGeoVolume *vol);
void PrintAllVolumes();

TGeoManager *geom = NULL;

void view_geom(std::string filename) {
  TGeoManager::Import(filename.c_str());
  geom = gGeoManager;
  geom->SetTopVisible(true); // Is the top layer visible?
  geom->DefaultColors();
  geom->SetVisOption(true);
  geom->SetVisLevel(199); // How many levels do we draw of the geometry

  // The whol world
  //geom->GetVolume("volWorld")->Draw("ogl");
  //geom->GetVolume("volWorld")->SetInvisible();
  // And the rock box
  //geom->GetVolume("rockBox_lv")->Draw("ogl");
  //geom->GetVolume("rockBox_lv")->SetInvisible();

  geom->GetVolume("volDetEnclosure")->SetLineColor(kBlack); // Make the hall visible
  geom->GetVolume("volDetEnclosure")->SetVisibility(true);
  Int_t trans = 50;

  // Some of the TMS volumes, useful for geometry inspection
  //geom->GetVolume("scinBoxlvTMS")->Draw("ogl");
  //geom->GetVolume("ModuleBoxvol")->Draw("ogl");
  //geom->GetVolume("modulelayervol1")->Draw("ogl");
  //geom->GetVolume("modulelayervol2")->Draw("ogl");
  //geom->GetVolume("modulelayervol3")->Draw("ogl");
  //geom->GetVolume("thickvol2TMS")->Draw("ogl");
  //geom->GetVolume("thickvolTMS")->Draw("ogl");
  //geom->GetVolume("thicklayervol")->Draw("ogl");
  //geom->GetVolume("thinvol2TMS")->Draw("ogl");
  //geom->GetVolume("thinvolTMS")->Draw("ogl");
  //geom->GetVolume("thinlayervol")->Draw("ogl");
  //geom->GetVolume("volTMS")->Draw("ogl");

  // Let's see if we can make everything visible

  //RecursiveVisible(geom->GetVolume("volTMS"));
  //RecursiveVisible(geom->GetVolume("volWorld"));

  RecursiveTransparency(geom->GetVolume("volWorld"), 10);
  RecursiveTransparency(geom->GetVolume("volTMS"), 80);
  //RecursiveTransparency(geom->GetVolume("volArgonCubeDetector_0"), 80);
  //geom->GetVolume("volDetEnclosure")->SetTransparency(80);

  // These will only exist in SAND geometries
  //RecursiveTransparency(geom->GetVolume("KLOEYokeBarrel_volume"), trans);

  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapAL_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapAR_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapBL_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapBR_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapCL_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapCR_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapDL_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOEYokeEndcapDR_volume"), trans);

  //RecursiveTransparency(geom->GetVolume("KLOESolenoidCryostatEndcapL_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOESolenoidCryostatEndcapR_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOESolenoidCryostatInnerWall_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOESolenoidCryostatOuterWall_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOESolenoidCoilShell_volume"), trans);
  //RecursiveTransparency(geom->GetVolume("KLOESolenoidCoil_volume"), trans);

  // Loop over the passive slabs
  /*TGeoVolume *ecalvol = geom->GetVolume("volECALPassiveSlab_0");
  int volno = 1;
  while (ecalvol != NULL) {
    TString temp = Form("volECALPassiveSlab_%i", volno);
    ecalvol = geom->GetVolume(temp.Data());
    if (ecalvol != NULL) ecalvol->SetTransparency(trans);
    volno++;
  }

  // Loop over the passive slabs
  ecalvol = geom->GetVolume("volECALActiveSlab_0");
  volno = 1;
  while (ecalvol != NULL) {
    TString temp = Form("volECALActiveSlab_%i", volno);
    ecalvol = geom->GetVolume(temp.Data());
    if (ecalvol != NULL) ecalvol->SetTransparency(trans);
    volno++;
  }*/


  //geom->GetVolume("ECAL_lv")->SetLineColor(kRed);
  //geom->GetVolume("ECAL_end_lv")->SetLineColor(kRed);
  //geom->GetVolume("kloe_calo_volume")->Draw("ogl");
  //geom->GetVolume("volKLOE")->Draw("ogl");
  //volMainDet_3DST
    //fullDetBox

  //geom->GetVolume("KLOEYokeEndcapAL_volume")->Draw("ogl");

  geom->GetVolume("volWorld")->Draw("ogl");

  // The support structure vlumes
  /*geom->GetVolume("elevatorBlock_lv")->SetLineColor(kBlack);
  geom->GetVolume("elevatorBlock_lv")->SetTransparency(trans);
  geom->GetVolume("craneRail1_lv")->SetLineColor(kBlack);
  geom->GetVolume("craneRail2_lv")->SetLineColor(kBlack);
  geom->GetVolume("egressHallway_lv")->SetLineColor(kBlack);
  geom->GetVolume("cryoTube1_lv")->SetLineColor(kBlack);
  geom->GetVolume("cryoTube2_lv")->SetLineColor(kBlack);
  geom->GetVolume("cryoTube3_lv")->SetLineColor(kBlack);
  geom->GetVolume("cryoTube4_lv")->SetLineColor(kBlack);
  geom->GetVolume("cryoTube5_lv")->SetLineColor(kBlack);*/

  // Example to not draw the support volume structures
  /*
  geom->GetVolume("elevatorBlock_lv")->SetVisibility(false);
  geom->GetVolume("craneRail1_lv")->SetVisibility(false);
  geom->GetVolume("craneRail2_lv")->SetVisibility(false);
  geom->GetVolume("egressHallway_lv")->SetVisibility(false);
  geom->GetVolume("cryoTube1_lv")->SetVisibility(false);
  geom->GetVolume("cryoTube2_lv")->SetVisibility(false);
  geom->GetVolume("cryoTube3_lv")->SetVisibility(false);
  geom->GetVolume("cryoTube4_lv")->SetVisibility(false);
  geom->GetVolume("cryoTube5_lv")->SetVisibility(false);
  geom->GetVolume("volCryostatSupportStructure")->SetVisibility(false);
  */

  //geom->GetVolume("sideplate_lv")->SetTransparency(90);
  //geom->GetVolume("endplate_lv")->SetTransparency(90);
  /*geom->GetVolume("sideplate_lv")->SetLineColor(kBlack);
  geom->GetVolume("endplate_lv")->SetLineColor(kBlack);
  geom->GetVolume("volTopInsulationR")->SetLineColor(kBlack);
  geom->GetVolume("volTopInsulationL")->SetLineColor(kBlack);
  geom->GetVolume("sideplate_lv")->SetTransparency(90);
  geom->GetVolume("endplate_lv")->SetTransparency(90);
  geom->GetVolume("volTopInsulationR")->SetTransparency(90);
  geom->GetVolume("volTopInsulationL")->SetTransparency(90);*/

  //geom->GetVolume("volDetEnclosure")->Draw("ogl");

  //geom->GetVolume("volCompositeWindow")->SetLineColor(kRed);
  //geom->GetVolume("volLArBath")->SetLineColor(kRed);
  //geom->GetVolume("volfullCryoTubeBox")->SetLineColor(kBlack);

  // The steel volumes
  // outer thick
  geom->GetVolume("thickvolTMS")->SetLineColor(kGreen);
  // middle thick
  geom->GetVolume("thickvol2TMS")->SetLineColor(kBlue);
  // outer thin
  geom->GetVolume("thinvolTMS")->SetLineColor(kYellow);
  // midle thin
  geom->GetVolume("thinvol2TMS")->SetLineColor(kMagenta);
  // outer double thick
  geom->GetVolume("doublevolTMS")->SetLineColor(kBlack);
  // middle double thick
  geom->GetVolume("doubelvol2TMS")->SetLineColor(kBlack);
  // The scintillator bar volume
  geom->GetVolume("scinBoxlvTMS")->SetLineColor(kRed);

  geom->GetVolume("modulelayervol1")->SetLineColor(kGreen);
  geom->GetVolume("modulelayervol2")->SetLineColor(kRed);
  geom->GetVolume("modulelayervol3")->SetLineColor(KBlue);
  //geom->GetVolume("thicklayervol")->SetLineColor(kRed);
  //geom->GetVolume("thinlayervol")->SetLineColor(kBlue);

  //geom->GetVolume("volfullCryoTubeBox")->SetVisContainers(true);
  //geom->GetVolume("volfullCryoTubeBox")->VisibleDaughters(true);
  //geom->GetVolume("volfullCryoTubeBox")->SetVisibility(true);

  //geom->GetVolume("volTMS")->Draw("ogl");
  //geom->GetVolume("thinvolTMS")->Draw("ogl");

  // Let's try do draw a track
  double beam_angle = 0.101;
  double p = 2;
  double px = p*sin(beam_angle);
  double py = 0;
  double pz = p*cos(beam_angle);
  double prodx = 0;
  double prody = 0;
  double prodz = 0;
  double time = 0;
  int pdg = 13;
  int status = 1;
  double E = sqrt(px*px+py*py+pz*pz+0.105*0.105);
  TParticle *part = new TParticle(pdg, status, 0, 0, 0, 0, px, py, pz, E, prodx, prody, prodz, time);
  int trackindex = geom->AddTrack(0, pdg, part);

  TVirtualGeoTrack *track = geom->GetTrack(trackindex);
  for (int i = 0; i < 1000; ++i) {
    track->AddPoint(0+i*5, 0+i*5, 0+i*100, 0);
  }
  track->SetLineColor(kRed);
  track->SetLineWidth(2);
  track->Print();
  track->SetLineColor(kRed);

  //PrintAllVolumes();

  //PrintAllVolumes("rockBox_lv");

  //TGeoVolume *vol = geom->GetVolume("volWorld");
  //vol->Draw("ogl");
  //PrintRecursive(geom->GetVolume("volWorld"));

  //PrintRecursive(geom->GetVolume("rockBox_lv"));

  //geom->GetVolume("volfullCryoTubeBox_PV")->Draw("ogl");
  //geom->GetVolume("volTMS")->Draw("ogl");
  //geom->GetVolume("scinBoxlvTMS_PV")->Draw("ogl");
  //geom->GetVolume("rockBox_lv_PV")->Draw("ogl");
  //geom->GetVolume("volWorld_PV")->Draw("ogl");
  //geom->GetVolume("volWorld")->Draw("ogl");
  //geom->GetVolume("rockBox_lv")->Draw("ogl");

  //TGLSAViewer *glsa = (TGLSAViewer *)gPad->GetViewer3D();
  // Draw each volume
  //geom->GetListOfVolumes()->Print();
  //geom->GetVolume("volCompositeWindow")->Draw("ogl");
  //geom->GetVolume("volTMS")->Draw("ogl");
  //geom->GetVolume("volWorld")->Draw("ogl");
  //TEveLine *line = new TEveLine;
  //for (int i = 0; i < 160; ++i) {
    //line->SetNextPoint(120*sin(0.2*i), 20*cos(0.2*i), 80-i);
  //}
  //gEve->AddElement(line);
  //track->Draw();

  //TObjArray* list = geom->GetVolume("volTMS")->GetListOfVolumes();
  //for (int i = 0; i < list->GetEntries(); ++i) {
    //TGeoVolume *vol = (TGeoVolume*)list->At(i);
    //vol->Print();
    ////vol->Draw("ogl");
  //}

  //geom->Draw("ogl");
  //geom->DrawTracks();
  //gPad->Modified();
  //gPad->Update();
  //track->Draw();

  //TGLViewer *v = (TGLViewer*)gPad->GetViewer3D();
  //v->SetStyle(TGLRnrCtx::kWireFrame);
  //v->SetSmoothPoints(kTRUE);
  //v->DrawGuides();
  //v->UpdateScene();
}

void PrintRecursive(TGeoVolume *vol) {
  if (vol == NULL) return;
  Int_t nd = vol->GetNdaughters();
  std::cout << vol->GetName() << ":" << std::endl;
  vol->Print();
  for (Int_t i=0; i<nd; i++) {
    std::cout << "Node " << i << ": " << vol->GetNode(i)->GetName() << " " << vol->GetNode(i)->GetVolume()->GetName() << std::endl;
    vol->GetNode(i)->GetVolume()->Print();
  }
}

void PrintAllVolumes() {
  TObjArray* list = geom->GetListOfVolumes();
  for (int i = 0; i < list->GetEntries(); ++i) {
    TGeoVolume *vol = (TGeoVolume*)list->At(i);
    vol->Print();
  }
}

void RecursiveInvisible(TGeoVolume *vol) {
  if (vol == NULL) return;
  vol->InvisibleAll();
  Int_t nd = vol->GetNdaughters();
  for (Int_t i=0; i<nd; i++) {
    RecursiveInvisible(vol->GetNode(i)->GetVolume());
  }
}

void RecursiveVisible(TGeoVolume *vol) {
  if (vol == NULL) return;
  vol->InvisibleAll(false);
  Int_t nd = vol->GetNdaughters();
  for (Int_t i=0; i<nd; i++) {
    RecursiveVisible(vol->GetNode(i)->GetVolume());
  }
}

void RecursiveTransparency(TGeoVolume *vol, Int_t transp) {
  if (vol == NULL) return;
  vol->SetTransparency(transp);
  vol->SetVisibility(true);
  Int_t nd = vol->GetNdaughters();
  for (Int_t i=0; i<nd; i++) {
    RecursiveTransparency(vol->GetNode(i)->GetVolume(), transp);
  }
}

void view_geom() {
  // Find these geometries in the ND_Production repo on github
  std::string filename = "nd_hall_with_lar_tms_nosand.gdml";
  view_geom(filename);
}

