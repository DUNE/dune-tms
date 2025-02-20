void validate_points(std::string filename) {
  TFile *f = new TFile(filename.c_str());
  TTree *t = (TTree*)f->Get("tree");
  gStyle->SetPalette(55);

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root","");
  canvname += "_validate.pdf";
  canv->Print(canvname+"[");

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","","colz");
  canv->Print(canvname);

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","muonReco==0","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","muonReco==0","colz");
  canv->Print(canvname);

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","muonReco==1","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","muonReco==1","colz");
  canv->Print(canvname);

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","muonReco==2","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","muonReco==2","colz");
  canv->Print(canvname);

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","muonReco==0 && muonStart == 1","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","muonReco==0 && muonStart == 1","colz");
  canv->Print(canvname);

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","muonReco==1 && muonStart == 1","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","muonReco==1 && muonStart == 1","colz");
  canv->Print(canvname);

  t->Draw("vtx[0]:vtx[2]>>h1(1000,0,2000,1000,-600,600)","muonReco==2 && muonStart == 1","colz");
  canv->Print(canvname);
  t->Draw("vtx[1]:vtx[2]>>h2(1000,0,2000,1000,-500,600)","muonReco==2 && muonStart == 1","colz");
  canv->Print(canvname);

  t->Draw("GlobalDeath[0]:GlobalDeath[2]>>h2(1000,0,3000,1000,-3000,1000)","","colz");
  canv->Print(canvname);
  t->Draw("GlobalDeath[1]:GlobalDeath[2]>>h2(1000,-100,3000,1000,-1000,1700)","","colz");
  canv->Print(canvname);

  t->Draw("GlobalDeath[0]:GlobalDeath[2]>>h2(1000,0,3000,1000,-3000,1000)","muonStart==0","colz");
  canv->Print(canvname);
  t->Draw("GlobalDeath[1]:GlobalDeath[2]>>h2(1000,-100,3000,1000,-1000,1700)","muonStart==0","colz");
  canv->Print(canvname);

  t->Draw("GlobalDeath[0]:GlobalDeath[2]>>h2(1000,0,3000,1000,-3000,1000)","muonStart==1","colz");
  canv->Print(canvname);
  t->Draw("GlobalDeath[1]:GlobalDeath[2]>>h2(1000,-100,3000,1000,-1000,1700)","muonStart==1","colz");
  canv->Print(canvname);

  t->Draw("GlobalDeath[0]:GlobalDeath[2]>>h2(1000,0,3000,1000,-3000,1000)","muonStart==2","colz");
  canv->Print(canvname);
  t->Draw("GlobalDeath[1]:GlobalDeath[2]>>h2(1000,-100,3000,1000,-1000,1700)","muonStart==2","colz");
  canv->Print(canvname);

  t->Draw("TMSDeath[0]:TMSDeath[2]>>h2(250,1000,2000,400,-400,400)","","colz");
  canv->Print(canvname);
  t->Draw("TMSDeath[1]:TMSDeath[2]>>h2(250,1000,2000,350,-300,50)","","colz");
  canv->Print(canvname);

  t->Draw("TMSDeath[0]:TMSDeath[2]>>h2(250,1000,2000,400,-400,400)","muonStart==0","colz");
  canv->Print(canvname);
  t->Draw("TMSDeath[1]:TMSDeath[2]>>h2(250,1000,2000,350,-300,50)","muonStart==0","colz");
  canv->Print(canvname);

  t->Draw("TMSDeath[0]:TMSDeath[2]>>h2(250,1000,2000,400,-400,400)","muonStart==1","colz");
  canv->Print(canvname);
  t->Draw("TMSDeath[1]:TMSDeath[2]>>h2(250,1000,2000,350,-300,50)","muonStart==1","colz");
  canv->Print(canvname);

  t->Draw("TMSDeath[0]:TMSDeath[2]>>h2(250,1000,2000,400,-400,400)","muonStart==2","colz");
  canv->Print(canvname);
  t->Draw("TMSDeath[1]:TMSDeath[2]>>h2(250,1000,2000,350,-300,50)","muonStart==2","colz");
  canv->Print(canvname);

  t->Draw("TMSBirth[0]:TMSBirth[2]>>h2(250,1000,2000,400,-400,400)","","colz");
  canv->Print(canvname);
  t->Draw("TMSBirth[1]:TMSBirth[2]>>h2(250,1000,2000,350,-300,50)","","colz");
  canv->Print(canvname);

  t->Draw("TMSBirth[0]:TMSBirth[2]>>h2(250,1000,2000,400,-400,400)","muonStart==0","colz");
  canv->Print(canvname);
  t->Draw("TMSBirth[1]:TMSBirth[2]>>h2(250,1000,2000,350,-300,50)","muonStart==0","colz");
  canv->Print(canvname);

  t->Draw("TMSBirth[0]:TMSBirth[2]>>h2(250,1000,2000,400,-400,400)","muonStart==1","colz");
  canv->Print(canvname);
  t->Draw("TMSBirth[1]:TMSBirth[2]>>h2(250,1000,2000,350,-300,50)","muonStart==1","colz");
  canv->Print(canvname);

  t->Draw("TMSBirth[0]:TMSBirth[2]>>h2(250,1000,2000,400,-400,400)","muonStart==2","colz");
  canv->Print(canvname);
  t->Draw("TMSBirth[1]:TMSBirth[2]>>h2(250,1000,2000,350,-300,50)","muonStart==2","colz");
  canv->Print(canvname);
  

  canv->Print(canvname+"]");
}
