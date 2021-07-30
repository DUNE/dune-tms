// Copied over from $ROOTSYS/tutorials/geom/geomDO.C

void RecursiveInvisible(TGeoVolume *vol);
void RecursiveVisible(TGeoVolume *vol);
void RecursiveTransparency(TGeoVolume *vol, Int_t transp);

void geom() {
  TFile *file = new TFile("geometry.root");
  TGeoManager *geom = (TGeoManager*)file->Get("geom");
  geom->SetTopVisible(true);
  geom->DefaultColors();
  geom->SetVisOption(0);
  geom->SetVisLevel(199);

  //geom->GetVolume("volWorld_PV")->Draw("ogl");
  //geom->GetVolume("scinBoxlvTMS_PV")->Draw("ogl");
  //geom->GetVolume("ModuleBoxvol_PV")->Draw("ogl");
  //geom->GetVolume("modulelayervol_PV")->Draw("ogl");
  //geom->GetVolume("thickvol2TMS_PV")->Draw("ogl");
  //geom->GetVolume("thickvolTMS_PV")->Draw("ogl");
  //geom->GetVolume("thicklayervol_PV")->Draw("ogl");
  //geom->GetVolume("thinvol2TMS_PV")->Draw("ogl");
  //geom->GetVolume("thinvolTMS_PV")->Draw("ogl");
  //geom->GetVolume("thinlayervol_PV")->Draw("ogl");

  // Let's see if we can make everything visible
  geom->GetListOfVolumes()->Print();
  //RecursiveVisible(geom->GetVolume("volTMS_PV"));
  //RecursiveVisible(geom->GetVolume("volWorld_PV"));
  RecursiveTransparency(geom->GetVolume("volWorld_PV"), 50);

  // The steel volumes
  // outer thick
  //geom->GetVolume("thickvolTMS_PV")->SetLineColor(kGreen);
  // middle thick
  //geom->GetVolume("thickvol2TMS_PV")->SetLineColor(kBlack);
  // outer thin
  //geom->GetVolume("thinvolTMS_PV")->SetLineColor(kYellow+3);
  // midle thin
  //geom->GetVolume("thinvol2TMS_PV")->SetLineColor(kGray);
  // The scintillator bar volume
  //geom->GetVolume("scinBoxlvTMS_PV")->SetLineColor(kRed);

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

  geom->GetVolume("volfullCryoTubeBox_PV")->SetVisContainers(true);
  geom->GetVolume("volfullCryoTubeBox_PV")->VisibleDaughters(true);
  geom->GetVolume("volfullCryoTubeBox_PV")->SetVisibility(true);

  //geom->GetVolume("volfullCryoTubeBox_PV")->Draw("ogl");
  //geom->GetVolume("volTMS_PV")->Draw("ogl");
  //geom->GetVolume("scinBoxlvTMS_PV")->Draw("ogl");
  //geom->GetVolume("rockBox_lv_PV")->Draw("ogl");
  geom->GetVolume("volWorld_PV")->Draw("ogl");
  geom->SetCurrentTrack(trackindex);
  TGLViewer *v = (TGLViewer*)gPad->GetViewer3D();
  //v->SetStyle(TGLRnrCtx::kWireFrame);
}


void RecursiveInvisible(TGeoVolume *vol) {
  vol->InvisibleAll();
  Int_t nd = vol->GetNdaughters();
  for (Int_t i=0; i<nd; i++) {
    RecursiveInvisible(vol->GetNode(i)->GetVolume());
  }
}

void RecursiveVisible(TGeoVolume *vol) {
  vol->InvisibleAll(false);
  Int_t nd = vol->GetNdaughters();
  for (Int_t i=0; i<nd; i++) {
    RecursiveVisible(vol->GetNode(i)->GetVolume());
  }
}

void RecursiveTransparency(TGeoVolume *vol, Int_t transp) {
  vol->SetTransparency(transp);
  vol->SetVisibility(true);
  Int_t nd = vol->GetNdaughters();
  for (Int_t i=0; i<nd; i++) {
    RecursiveTransparency(vol->GetNode(i)->GetVolume(), transp);
  }
}

/*
 == Volume: scinBoxlvTMS_PV type TGeoVolume positioned 48 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =    17.71000
    dY =  1548.00000
    dZ =     5.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Scintillator    Aeff=11.0865 Zeff=5.58 rho=6.55359e+18 radlen=417.968 intlen=688.732 index=4
   Element #0 : C  Z=  6.00 A= 12.00 w= 0.906
   Element #1 : C  Z=  6.00 A= 13.00 w= 0.010
   Element #2 : H  Z=  1.00 A=  1.01 w= 0.084
   Element #3 : H  Z=  1.00 A=  2.01 w= 0.000


   == Volume: ModuleBoxvol_PV type TGeoVolume positioned 56 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =   850.08000
    dY =  1548.00000
    dZ =     5.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Air    Aeff=14.666 Zeff=7.31201 rho=7.64585e+15 radlen=301550 intlen=697784 index=2
   Element #0 : N  Z=  7.00 A= 14.00 w= 0.778
   Element #1 : N  Z=  7.00 A= 15.00 w= 0.003
   Element #2 : O  Z=  8.00 A= 15.99 w= 0.209
   Element #3 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #4 : O  Z=  8.00 A= 18.00 w= 0.000
   Element #5 : AR  Z= 18.00 A= 35.97 w= 0.000
   Element #6 : AR  Z= 18.00 A= 37.96 w= 0.000
   Element #7 : AR  Z= 18.00 A= 39.96 w= 0.009

    == Volume: modulelayervol_PV type TGeoVolume positioned 108 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =  3518.00000
    dY =  2511.00000
    dZ =    20.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Air    Aeff=14.666 Zeff=7.31201 rho=7.64585e+15 radlen=301550 intlen=697784 index=2
   Element #0 : N  Z=  7.00 A= 14.00 w= 0.778
   Element #1 : N  Z=  7.00 A= 15.00 w= 0.003
   Element #2 : O  Z=  8.00 A= 15.99 w= 0.209
   Element #3 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #4 : O  Z=  8.00 A= 18.00 w= 0.000
   Element #5 : AR  Z= 18.00 A= 35.97 w= 0.000
   Element #6 : AR  Z= 18.00 A= 37.96 w= 0.000
   Element #7 : AR  Z= 18.00 A= 39.96 w= 0.009

    == Volume: thickvol2TMS_PV type TGeoVolume positioned 1 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =  1749.00000
    dY =  2511.00000
    dZ =    20.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture SteelTMS    Aeff=55.7027 Zeff=25.938 rho=4.89959e+19 radlen=17.6625 intlen=170.18 index=3
   Element #0 : FE  Z= 26.00 A= 53.94 w= 0.058
   Element #1 : FE  Z= 26.00 A= 55.93 w= 0.913
   Element #2 : FE  Z= 26.00 A= 56.94 w= 0.021
   Element #3 : FE  Z= 26.00 A= 57.93 w= 0.003
   Element #4 : SI  Z= 14.00 A= 27.98 w= 0.004
   Element #5 : SI  Z= 14.00 A= 28.98 w= 0.000
   Element #6 : SI  Z= 14.00 A= 29.97 w= 0.000
   Element #7 : C  Z=  6.00 A= 12.00 w= 0.000
   Element #8 : C  Z=  6.00 A= 13.00 w= 0.000

    == Volume: thickvolTMS_PV type TGeoVolume positioned 2 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =   874.50000
    dY =  2511.00000
    dZ =    20.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture SteelTMS    Aeff=55.7027 Zeff=25.938 rho=4.89959e+19 radlen=17.6625 intlen=170.18 index=3
   Element #0 : FE  Z= 26.00 A= 53.94 w= 0.058
   Element #1 : FE  Z= 26.00 A= 55.93 w= 0.913
   Element #2 : FE  Z= 26.00 A= 56.94 w= 0.021
   Element #3 : FE  Z= 26.00 A= 57.93 w= 0.003
   Element #4 : SI  Z= 14.00 A= 27.98 w= 0.004
   Element #5 : SI  Z= 14.00 A= 28.98 w= 0.000
   Element #6 : SI  Z= 14.00 A= 29.97 w= 0.000
   Element #7 : C  Z=  6.00 A= 12.00 w= 0.000
   Element #8 : C  Z=  6.00 A= 13.00 w= 0.000

    == Volume: thicklayervol_PV type TGeoVolume positioned 63 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =  3518.00000
    dY =  2511.00000
    dZ =    20.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Air    Aeff=14.666 Zeff=7.31201 rho=7.64585e+15 radlen=301550 intlen=697784 index=2
   Element #0 : N  Z=  7.00 A= 14.00 w= 0.778
   Element #1 : N  Z=  7.00 A= 15.00 w= 0.003
   Element #2 : O  Z=  8.00 A= 15.99 w= 0.209
   Element #3 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #4 : O  Z=  8.00 A= 18.00 w= 0.000
   Element #5 : AR  Z= 18.00 A= 35.97 w= 0.000
   Element #6 : AR  Z= 18.00 A= 37.96 w= 0.000
   Element #7 : AR  Z= 18.00 A= 39.96 w= 0.009

    == Volume: thinvol2TMS_PV type TGeoVolume positioned 1 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =  1749.00000
    dY =  2511.00000
    dZ =     7.50000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture SteelTMS    Aeff=55.7027 Zeff=25.938 rho=4.89959e+19 radlen=17.6625 intlen=170.18 index=3
   Element #0 : FE  Z= 26.00 A= 53.94 w= 0.058
   Element #1 : FE  Z= 26.00 A= 55.93 w= 0.913
   Element #2 : FE  Z= 26.00 A= 56.94 w= 0.021
   Element #3 : FE  Z= 26.00 A= 57.93 w= 0.003
   Element #4 : SI  Z= 14.00 A= 27.98 w= 0.004
   Element #5 : SI  Z= 14.00 A= 28.98 w= 0.000
   Element #6 : SI  Z= 14.00 A= 29.97 w= 0.000
   Element #7 : C  Z=  6.00 A= 12.00 w= 0.000
   Element #8 : C  Z=  6.00 A= 13.00 w= 0.000

    == Volume: thinvolTMS_PV type TGeoVolume positioned 2 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =   874.50000
    dY =  2511.00000
    dZ =     7.50000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture SteelTMS    Aeff=55.7027 Zeff=25.938 rho=4.89959e+19 radlen=17.6625 intlen=170.18 index=3
   Element #0 : FE  Z= 26.00 A= 53.94 w= 0.058
   Element #1 : FE  Z= 26.00 A= 55.93 w= 0.913
   Element #2 : FE  Z= 26.00 A= 56.94 w= 0.021
   Element #3 : FE  Z= 26.00 A= 57.93 w= 0.003
   Element #4 : SI  Z= 14.00 A= 27.98 w= 0.004
   Element #5 : SI  Z= 14.00 A= 28.98 w= 0.000
   Element #6 : SI  Z= 14.00 A= 29.97 w= 0.000
   Element #7 : C  Z=  6.00 A= 12.00 w= 0.000
   Element #8 : C  Z=  6.00 A= 13.00 w= 0.000

    == Volume: thinlayervol_PV type TGeoVolume positioned 42 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =  3518.00000
    dY =  2511.00000
    dZ =     7.50000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Air    Aeff=14.666 Zeff=7.31201 rho=7.64585e+15 radlen=301550 intlen=697784 index=2
   Element #0 : N  Z=  7.00 A= 14.00 w= 0.778
   Element #1 : N  Z=  7.00 A= 15.00 w= 0.003
   Element #2 : O  Z=  8.00 A= 15.99 w= 0.209
   Element #3 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #4 : O  Z=  8.00 A= 18.00 w= 0.000
   Element #5 : AR  Z= 18.00 A= 35.97 w= 0.000
   Element #6 : AR  Z= 18.00 A= 37.96 w= 0.000
   Element #7 : AR  Z= 18.00 A= 39.96 w= 0.009

    == Volume: volTMS_PV type TGeoVolume positioned 200 times
*** Shape TGeoBBox: TGeoBBox ***
    dX =  3518.00000
    dY =  3450.00000
    dZ =  3525.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Air    Aeff=14.666 Zeff=7.31201 rho=7.64585e+15 radlen=301550 intlen=697784 index=2
   Element #0 : N  Z=  7.00 A= 14.00 w= 0.778
   Element #1 : N  Z=  7.00 A= 15.00 w= 0.003
   Element #2 : O  Z=  8.00 A= 15.99 w= 0.209
   Element #3 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #4 : O  Z=  8.00 A= 18.00 w= 0.000
   Element #5 : AR  Z= 18.00 A= 35.97 w= 0.000
   Element #6 : AR  Z= 18.00 A= 37.96 w= 0.000
   Element #7 : AR  Z= 18.00 A= 39.96 w= 0.009

    == Volume: volDetEnclosure_PV type TGeoVolume positioned 8 times
*** TGeoCompositeShape : name =
 Bounding box:
*** Shape name: TGeoBBox ***
    dX = 27587.50000
    dY = 12215.00000
    dZ = 14927.40000
    origin: x=-2892.50000 y= 4595.00000 z= 5022.60000
Mixture Air    Aeff=14.666 Zeff=7.31201 rho=7.64585e+15 radlen=301550 intlen=697784 index=2
   Element #0 : N  Z=  7.00 A= 14.00 w= 0.778
   Element #1 : N  Z=  7.00 A= 15.00 w= 0.003
   Element #2 : O  Z=  8.00 A= 15.99 w= 0.209
   Element #3 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #4 : O  Z=  8.00 A= 18.00 w= 0.000
   Element #5 : AR  Z= 18.00 A= 35.97 w= 0.000
   Element #6 : AR  Z= 18.00 A= 37.96 w= 0.000
   Element #7 : AR  Z= 18.00 A= 39.96 w= 0.009

    == Volume: volWorld_PV type TGeoVolume positioned 1 times
*** Shape TGeoBBox: TGeoBBox ***
    dX = 300000.00000
    dY = 300000.00000
    dZ = 300000.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture NoGas    Aeff=39.9477 Zeff=18 rho=6.24151e-07 radlen=1.9549e+27 intlen=1.19646e+28 index=0
   Element #0 : AR  Z= 18.00 A= 35.97 w= 0.003
   Element #1 : AR  Z= 18.00 A= 37.96 w= 0.001
   Element #2 : AR  Z= 18.00 A= 39.96 w= 0.996
 == Volume: rockBox_lv_PV type TGeoVolume positioned 2 times
*** TGeoCompositeShape : name =
 Bounding box:
*** Shape name: TGeoBBox ***
    dX = 84585.00000
    dY = 67820.00000
    dZ = 139775.00000
    origin: x=    0.00000 y=    0.00000 z=    0.00000
Mixture Rock    Aeff=24.2904 Zeff=11.9245 rho=1.76011e+19 radlen=90.1394 intlen=345.029 index=1
   Element #0 : SI  Z= 14.00 A= 27.98 w= 0.227
   Element #1 : SI  Z= 14.00 A= 28.98 w= 0.012
   Element #2 : SI  Z= 14.00 A= 29.97 w= 0.008
   Element #3 : O  Z=  8.00 A= 15.99 w= 0.493
   Element #4 : O  Z=  8.00 A= 17.00 w= 0.000
   Element #5 : O  Z=  8.00 A= 18.00 w= 0.001
   Element #6 : FE  Z= 26.00 A= 53.94 w= 0.005
   Element #7 : FE  Z= 26.00 A= 55.93 w= 0.084
   Element #8 : FE  Z= 26.00 A= 56.94 w= 0.002
   Element #9 : FE  Z= 26.00 A= 57.93 w= 0.000
   Element #10 : AL  Z= 13.00 A= 26.98 w= 0.054
   Element #11 : MG  Z= 12.00 A= 23.98 w= 0.023
   Element #12 : MG  Z= 12.00 A= 24.99 w= 0.003
   Element #13 : MG  Z= 12.00 A= 25.98 w= 0.003
   Element #14 : C  Z=  6.00 A= 12.00 w= 0.035
   Element #15 : C  Z=  6.00 A= 13.00 w= 0.000
   Element #16 : CA  Z= 20.00 A= 39.96 w= 0.026
   Element #17 : CA  Z= 20.00 A= 41.96 w= 0.000
   Element #18 : CA  Z= 20.00 A= 42.96 w= 0.000
   Element #19 : CA  Z= 20.00 A= 43.96 w= 0.001
   Element #20 : CA  Z= 20.00 A= 45.95 w= 0.000
   Element #21 : CA  Z= 20.00 A= 47.95 w= 0.000
   Element #22 : S  Z= 16.00 A= 31.97 w= 0.018
   Element #23 : S  Z= 16.00 A= 32.97 w= 0.000
   Element #24 : S  Z= 16.00 A= 33.97 w= 0.001
   Element #25 : S  Z= 16.00 A= 35.97 w= 0.000
   Element #26 : NA  Z= 11.00 A= 22.99 w= 0.004
   Element #27 : P  Z= 15.00 A= 30.97 w= 0.000
*/
