// python code was very slow...
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"

#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4PrimaryVertex.h"

const double mass_mu = 105.6583755; // Muon mass in MeV/c2
const double extra_trk_length = 24.35;
const double LAr_density = 1.3954;
const double TMS_Scint_Width = 4; // TMS scintillator width (4 cm)
const double TMS_Thin_Layer_z_end = 949; // TMS region that is thin iron layer -> 5.5cm gap
const double TMS_Transition_Layer_z_end = 949.3; // TMS region that is thin iron layer -> 5.5cm gap
const double TMS_Thin_Layer_gap = 5.5; // Gap for TMS region that is thin iron layer
const double TMS_Transition_Layer_gap = 9.5; // Gap for TMS region that is between thin and thick regions
const double TMS_Thick_Layer_gap = 8.0; // Gap for TMS region that is thick iron layer
const int TMS_Max_Planes = 8; // Number of planes in TMS

bool dumpSSRITree(std::string filename, std::string output_filename) {
  std::cout << "Got " << filename << ", writing to " << output_filename << std::endl;

  // Offset offsets (all?) the vertex position in {x,y,z}, don't know why
  const double offset[] = { 0., 5.5, 411. };
  // Fiducial volume cuts? Never used
  const double fvLo[] = { -300., -100., 50. };
  const double fvHi[] = { 300., 100., 450. };
  // Never used
  const double collarLo[] = { -320., -120., 30. };
  const double collarHi[] = { 320., 120., 470. };

  const int nplanes_cut = 8;

  TFile *input = new TFile(filename.c_str());
  // The EDepSim events
  TTree *events = (TTree*)input->Get("EDepSimEvents");
  // The generator pass-through information
  TTree *gRoo = (TTree*)input->Get("DetSimPassThru/gRooTracker");
  // Get the detector geometry
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);
  // For the true event we just care about the neutrino
  // GENIE saves max 100 particles in gRooTracker format
  const int kNPmax = 100;
  int True_NuPDG[kNPmax];
  double True_NuE[kNPmax][4];
  gRoo->SetBranchAddress("StdHepPdg", True_NuPDG);
  gRoo->SetBranchAddress("StdHepP4", True_NuE);

  int N_entries = events->GetEntries();
  int ient = 0;

  //////////////////////////////////////////////////
  // Setup the output tree
  TFile *output = new TFile(output_filename.c_str(), "RECREATE" );
  TTree *tree_out = new TTree("tree", "tree");

  std::string ReactionString;
  int ievt;
  bool SegInLAr;
  bool SegInTMS;
  float E_nu;
  int PDG_nu;
  int LepPDG;
  float LepKE;
  float RMMS_KE;
  float LepE;
  float MuonExitKE;
  int MuonReco;
  float MuonScintLen;
  float MuonLArLen;
  float MuonScintLen_Old;
  float MuonScintEnergy;
  int nFinalState;
  float Lepton[3];
  float Vertex[3];
  float LeptonDeath[3];
  float MuonDeath[3];
  float MuonBirth[3];
  float LepRMMS[3];
  float MuonExitPt[3];
  float MuonExitMom[3];
  float MuonTMSEntryMom[3];
  const int MAX_PARTICLES = 100;
  int PDGFinalState[MAX_PARTICLES];
  float PxFinalState[MAX_PARTICLES];
  float PyFinalState[MAX_PARTICLES];
  float PzFinalState[MAX_PARTICLES];
  float EFinalState[MAX_PARTICLES];
  std::vector<float> xPos;
  std::vector<float> yPos;
  std::vector<float> zPos;
  std::vector<float> tPos;
  std::vector<float> xPosLep;
  std::vector<float> yPosLep;
  std::vector<float> zPosLep;
  std::vector<float> tPosLep;

  tree_out->Branch("reac", &ReactionString);
  tree_out->Branch("ievt", &ievt, "ievt/I");
  tree_out->Branch("Ev", &E_nu, "Ev/F");
  tree_out->Branch("PDGv", &PDG_nu, "PDGv/I");
  tree_out->Branch("p3lep", Lepton, "p3lep[3]/F");
  tree_out->Branch("vtx", Vertex, "vtx[3]/F");
  tree_out->Branch("lepDeath", LeptonDeath, "lepDeath[3]/F");
  tree_out->Branch("muonDeath", MuonDeath, "muonDeath[3]/F");
  tree_out->Branch("muonBirth", MuonBirth, "muonBirth[3]/F");
  tree_out->Branch("lepRMMS", LepRMMS, "lepRMMS[3]/F");
  tree_out->Branch("lepPdg", &LepPDG, "lepPDG/I");
  tree_out->Branch("lepKE", &LepKE, "lepKE/F");
  tree_out->Branch("rmmsKE", &RMMS_KE, "rmmsKE/F");
  tree_out->Branch("lepE", &LepE, "lepE/F");
  tree_out->Branch("muonExitPt", MuonExitPt, "muonExitPt[3]/F");
  tree_out->Branch("muonExitMom", MuonExitMom, "muonExitMom[3]/F");
  tree_out->Branch("muonExitKE", &MuonExitKE, "muonExitKE/F");
  tree_out->Branch("muonTMSEntryMom", &MuonTMSEntryMom, "muonTMSEntryMom[3]/F");
  tree_out->Branch("muonReco", &MuonReco, "muonReco/I");
  tree_out->Branch("muScintLen", &MuonScintLen, "muScintLen/F");
  tree_out->Branch("muLArLen", &MuonLArLen, "muLArLen/F");
  tree_out->Branch("muScintLenOld", &MuonScintLen_Old, "muScintLenOld/F");
  tree_out->Branch("muScintEnergy", &MuonScintEnergy, "muScintEnergy/F");
  tree_out->Branch("LArSeg", &SegInLAr, "LArSeg/O");
  tree_out->Branch("TMSSeg", &SegInTMS, "TMSSeg/O");
  tree_out->Branch("nFS", &nFinalState, "nFS/I");
  tree_out->Branch("fsPdg", PDGFinalState, "fsPdg[nFS]/I");
  tree_out->Branch("fsPx", PxFinalState, "fsPx[nFS]/F");
  tree_out->Branch("fsPy", PyFinalState, "fsPy[nFS]/F");
  tree_out->Branch("fsPz", PzFinalState, "fsPz[nFS]/F");
  tree_out->Branch("fsE", EFinalState, "fsE[nFS]/F");
  tree_out->Branch("xpt", &xPos);
  tree_out->Branch("ypt", &yPos);
  tree_out->Branch("zpt", &zPos);
  tree_out->Branch("tpt", &tPos);
  tree_out->Branch("xptLep", &xPosLep);
  tree_out->Branch("yptLep", &yPosLep);
  tree_out->Branch("zptLep", &zPosLep);
  tree_out->Branch("tptLep", &tPosLep);
  // Finished the output branches
  //////////////////////////////////////////////////

  // Now loop
  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();
  for (int i = 0; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    E_nu = True_NuE[0][3];
    PDG_nu = True_NuPDG[0];

    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
      std::cout << "  E_nu: " << E_nu << " PDG: " << PDG_nu << std::endl;
    }

    // Loop over the primary vertices
    for (TG4PrimaryVertexContainer::iterator it = event->Primaries.begin(); it != event->Primaries.end(); ++it) {
      // The vertex
      TG4PrimaryVertex vtx = *it;

      // Now reset the variables
      ievt = LepPDG = MuonReco = -999;
      LepKE = RMMS_KE = LepE = MuonExitKE = MuonScintLen = MuonLArLen = MuonScintLen_Old = MuonScintEnergy = -999;
      nFinalState = 0;
      for (int j = 0; j < 3; ++j) {
        Lepton[j] = Vertex[j] = LeptonDeath[j] = MuonDeath[j] = MuonBirth[j] = LepRMMS[j] = MuonExitPt[j] = MuonExitMom[j] = MuonTMSEntryMom[j] = -999;
      }
      for (int j = 0; j < MAX_PARTICLES; ++j) {
        PDGFinalState[j] = PxFinalState[j] = PyFinalState[j] = PzFinalState[j] = EFinalState[j] = -999;
      }
      xPos.clear();
      yPos.clear();
      zPos.clear();
      tPos.clear();
      xPosLep.clear();
      yPosLep.clear();
      zPosLep.clear();
      tPosLep.clear();
      ReactionString.clear();
      // Done resetting variables

      ReactionString = vtx.GetReaction();
      ievt = i;
      // Technically a TLorentzVector
      for (int j = 0; j < 3; ++j) Vertex[j] = (vtx.GetPosition()[j])/10. - offset[j]; // in cm

      // Find the lepton
      int LepIndex = -1;
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;
      // Loop over the particles in the vertex
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
        // Get some fundamentals about the particle
        int pdg = particle.GetPDGCode();

        PDGFinalState[nFinalState] = pdg;
        PxFinalState[nFinalState] = particle.GetMomentum().Px();
        PyFinalState[nFinalState] = particle.GetMomentum().Py();
        PzFinalState[nFinalState] = particle.GetMomentum().Pz();
        EFinalState[nFinalState] = particle.GetMomentum().E();

        // Find the lepton with matching PDG to incoming neutrino
        // If NC, PDG matches
        if (PDG_nu == pdg) LepIndex = nFinalState;
        // If neutrino, look for charged lepton PDG of PDG_nu-1
        else if (PDG_nu > 0 && PDG_nu == pdg+1) LepIndex = nFinalState;
        // If anti-neutrino, look for charged lepton PDG of PDG_nu+1
        else if (PDG_nu < 0 && PDG_nu == pdg-1) LepIndex = nFinalState;
        //if (abs(pdg) == 13) LepIndex = nFinalState;

        // Increment final state particle counter
        nFinalState++;
      }

      // Look only at muons for now
      if (abs(PDGFinalState[LepIndex]) != 13) continue;

      // Check that we had an outgoing lepton
      if (LepIndex == -1) {
        std::cerr << "No lepton found in event " << ievt << std::endl;
      }

      // Now extract the Lepton kinematics
      TG4PrimaryParticle lepton = particles[LepIndex];
      for (int j = 0; j < 3; ++j) Lepton[j] = lepton.GetMomentum()[j];
      LepPDG = lepton.GetPDGCode();
      LepKE = lepton.GetMomentum().E()-lepton.GetMomentum().M();
      LepE = EFinalState[LepIndex];

      // Also save the PrimaryId?
      int LepTrackId = particles[LepIndex].GetTrackId();
      //std::cout << "leptrackid: " << LepTrackId << std::endl;


      // If there is a muon, determine how to reconstruct its momentum and charge
      bool Exit = false;
      bool InTMS = false;

      // Get the trajectory of that lepton trajector
      TG4Trajectory LepTraj = event->Trajectories[LepIndex];
      // Loop over the points that the lepton trajectory left
      std::vector<TG4TrajectoryPoint> LepPoints = LepTraj.Points;
      TLorentzVector PreviousPoint(-999, -999, -999, -999);
      int index = 0;
      int total = LepPoints.size();
      //std::cout << "***" << std::endl;
      //std::cout << "Event " << ievt << " lepton points" << std::endl;
      for (TG4Trajectory::TrajectoryPoints::iterator jt = LepPoints.begin(); jt != LepPoints.end(); ++jt) {
        TG4TrajectoryPoint pt = *jt;
        TLorentzVector Point = pt.GetPosition();
        TGeoNode *node = geom->FindNode(Point.X(), Point.Y(), Point.Z());
        if (node == NULL) {
          std::cerr << "Couldn't find node for this point" << std::endl;
          continue;
        }
        std::string VolumeName = node->GetName();
        bool active = false;
        //std::cout << index << "/" << total << ", name: " << VolumeName << " (x,y,z)=(" << Point.X()/10.-offset[0] << ", " << Point.Y()/10.-offset[1] << ", " << Point.Z()/10.-offset[2] << ") dist in z: " << (Point.Z()/10.)-(PreviousPoint.Z()/10.) << " (px,py,pz)=(" << pt.GetMomentum().Px() << ", " << pt.GetMomentum().Py() << ", " << pt.GetMomentum().Pz() << ")" << std::endl;

        // Easiest to just do a string comparison, unfortunately
        // If in active LAr volume or pixels, update the exit points
        if (VolumeName.find("LAr") != std::string::npos || 
            VolumeName.find("PixelPlane") != std::string::npos || 
            VolumeName.find("sPlane") != std::string::npos) {
          // What momentum did it have when it exited the LAr
          MuonExitMom[0] = pt.GetMomentum().Px();
          MuonExitMom[1] = pt.GetMomentum().Py();
          MuonExitMom[2] = pt.GetMomentum().Pz();
          // Assume it's a muon to get a KE estimate
          MuonExitKE = sqrt(pt.GetMomentum().Mag2()+mass_mu*mass_mu) - mass_mu;
          MuonExitPt[0] = Point.X()/10.-offset[0];
          MuonExitPt[1] = Point.Y()/10.-offset[1];
          MuonExitPt[2] = Point.Z()/10.-offset[2];
        // If not in tracked LAr
        } else {

          // Check if it's inside the TMS
          // Read what the kinetic energy and momentum when it is when entering the TMS
          if (!InTMS && (VolumeName.find("RMMS") != std::string::npos || 
                         VolumeName.find("modulelayer") != std::string::npos)) {
            RMMS_KE = sqrt(pt.GetMomentum().Mag2()+mass_mu*mass_mu) - mass_mu;
            MuonTMSEntryMom[0] = pt.GetMomentum().Px();
            MuonTMSEntryMom[1] = pt.GetMomentum().Py();
            MuonTMSEntryMom[2] = pt.GetMomentum().Pz();
            InTMS = true;
          }
        }
        // Update the previous point
        PreviousPoint = Point;
        index++;
      }

      // Get the last point to see where the trajectory went
      TLorentzVector EndPoint = LepPoints.back().GetPosition();
      TGeoNode *EndNode = geom->FindNode(EndPoint.X(), EndPoint.Y(), EndPoint.Z());
      if (EndNode == NULL) {
        std::cerr << "Found no end node for this event" << std::endl;
      } else {
        LeptonDeath[0] = EndPoint.X()/10.-offset[0];
        LeptonDeath[1] = EndPoint.Y()/10.-offset[1];
        LeptonDeath[2] = EndPoint.Z()/10.-offset[2];
        std::string EndVolName = EndNode->GetName();
        if (EndVolName.find("LArActive") != std::string::npos) MuonReco = 1;
        else if (EndVolName.find("RMMS") != std::string::npos) MuonReco = 2;
        else MuonReco = 0;
      }

      // Remember if we've seen hits in the LAr or the TMS
      SegInLAr = false;
      SegInTMS = false;

      for (TG4HitSegmentDetectors::iterator jt = event->SegmentDetectors.begin(); jt != event->SegmentDetectors.end(); ++jt) {

        // Get the name of the active detector
        std::string DetString = (*jt).first;

        // Remember what this segment is in
        bool ThisSegInLAr = false;
        bool ThisSegInTMS = false;

        // If in ArgonCube
        if (DetString == "ArgonCube") {
          ThisSegInLAr = true;
        // If in RMMS
        } else if (DetString == "rmmsvol") {
          ThisSegInTMS = true;
        }

        // The container for all of the segments
        TG4HitSegmentContainer hitsegments = (*jt).second;

        // Save the lepton only hits
        TG4HitSegmentContainer LeptonHits;
        // Save the other hits
        TG4HitSegmentContainer AllHits;

        // Then loop over each of hit segments in this container
        // Now look at the individual hits in ArgonCube and TMS for this event
        //std::cout << "***" << std::endl;
        for (TG4HitSegmentContainer::iterator kt = hitsegments.begin(); kt != hitsegments.end(); ++kt) {
          // The current segment
          TG4HitSegment seg = (*kt);
          //std::cout << seg.GetPrimaryId() << std::endl;

          // Loop over the contributors to this track segment
          for (TG4HitSegment::Contributors::iterator ct = seg.Contrib.begin(); ct != seg.Contrib.end(); ++ct) {
            int Contributor = *ct;
            //int Contributor = 0;
            // Get the trajectory that produced this hit
            TG4Trajectory Trajectory = event->Trajectories[Contributor];
            int ParentId = Trajectory.GetParentId();
            int pdg = Trajectory.GetPDGCode();
            int TrackId = Trajectory.GetTrackId();
            // If from the fundamental interaction vertex and the lepton
            if (ParentId == -1 && Contributor == TrackId && TrackId == LepTrackId) {
              // Check if it's a lepton matching the neutrino PDG
              if ((PDG_nu > 0 && PDG_nu == pdg+1) || 
                   PDG_nu < 0 && PDG_nu == pdg-1) {
              //if (abs(pdg) == 13) {
                LeptonHits.push_back(seg);
              } else {
                AllHits.push_back(seg);
              }
              // Not from the primary interaction
            } else {
              AllHits.push_back(seg);
            }
          }
        }

          //std::cout << "***" << std::endl;
        // Now we have all the hits and all the lepton hits
        // Save them into the positions
        for (TG4HitSegmentContainer::iterator kt = LeptonHits.begin(); kt != LeptonHits.end(); ++kt) {
          TG4HitSegment seg = *kt;
          xPosLep.push_back(seg.GetStart()[0]/10.-offset[0]);
          yPosLep.push_back(seg.GetStart()[1]/10.-offset[1]);
          zPosLep.push_back(seg.GetStart()[2]/10.-offset[2]);
          tPosLep.push_back(seg.GetStart()[3]);
        }
        for (TG4HitSegmentContainer::iterator kt = AllHits.begin(); kt != AllHits.end(); ++kt) {
          TG4HitSegment seg = *kt;
          xPos.push_back(seg.GetStart()[0]/10.-offset[0]);
          yPos.push_back(seg.GetStart()[1]/10.-offset[1]);
          zPos.push_back(seg.GetStart()[2]/10.-offset[2]);
          tPos.push_back(seg.GetStart()[3]);
        }

        // If this segment is in the LAr
        if (ThisSegInLAr) {
          // Require at least three hit segments in the LAr for our lepton, else don't count it
          if (LeptonHits.size() < 3) {
            continue;
          }
          // Remember that we have a valid segment in ArgonCube
          SegInLAr = true;
          TVector3 LArStart(LeptonHits.front().GetStart()[0]/10.-offset[0],
              LeptonHits.front().GetStart()[1]/10.-offset[1],
              LeptonHits.front().GetStart()[2]/10.-offset[2]);
          TVector3 LArEnd(LeptonHits.back().GetStart()[0]/10.-offset[0],
              LeptonHits.back().GetStart()[1]/10.-offset[1],
              LeptonHits.back().GetStart()[2]/10.-offset[2]);
          double LAr_dist = (LArEnd-LArStart).Mag();
          double LAr_dist_z = LArEnd.z()-LArStart.z();
          // Track length
          double LAr_Track_Length = LAr_density*LAr_dist;
          MuonLArLen = LAr_Track_Length;

          //for (int dim = 0; dim < 3; ++dim) MuonExitPt[dim] = LeptonHits.back().GetStart()[dim]/10. - offset[dim];

          // If this segment is in the TMS
        } else if (ThisSegInTMS) {
          // Need more than one hit
          if (LeptonHits.size() < 2) {
            continue;
          }
          // Remember that we have a valid segment in the TMS
          SegInTMS = true;
          // Save the muon birth point in the TMS
          for (int dim = 0; dim < 3; ++dim) MuonBirth[dim] = LeptonHits.front().GetStart()[dim]/10. - offset[dim];
          // Now count up the deposited energy in the lepton track
          double TMS_EnergyDeposit = 0;
          double TMS_TrackLength = 0;
          // Keep the previous hit position
          TVector3 TMS_HitStart_Prev;
          TVector3 TMS_HitEnd_Prev;
          // How many planes did we traverse
          int TMS_total_planes = 0;
          for (TG4HitSegmentContainer::iterator kt = LeptonHits.begin(); kt != LeptonHits.end(); ++kt) {
            // Get the hit segment
            TG4HitSegment TMS_seg = *kt;
            //std::cout << "***" << std::endl;
            //std::cout << "Tracklength for hit in edepsim: " << TMS_seg.GetTrackLength() << std::endl;
            //std::cout << " TMS hit Lep (x,y,z) = ";
            //for (int dim = 0; dim < 3; ++ dim) std::cout << TMS_seg.GetStart()[dim]/10-offset[dim] << " to " << TMS_seg.GetStop()[dim]/10-offset[dim] << ", ";

            
            // Incremement the energy deposit
            TMS_EnergyDeposit += TMS_seg.GetEnergyDeposit();

            TVector3 TMS_HitStart(TMS_seg.GetStart()[0]/10.-offset[0],
                TMS_seg.GetStart()[1]/10.-offset[1],
                TMS_seg.GetStart()[2]/10.-offset[2]);
            TVector3 TMS_HitEnd(TMS_seg.GetStop()[0]/10.-offset[0],
                TMS_seg.GetStop()[1]/10.-offset[1],
                TMS_seg.GetStop()[2]/10.-offset[2]);

            // Subtract off the previous hit
            //std::cout << " z dist: " << TMS_HitStart.Z()-TMS_HitStart_Prev.Z() << std::endl;

            // Check that we've moved hit
            if (TMS_HitStart_Prev.Mag() == 0) {
              TMS_HitStart_Prev = TMS_HitStart;
              TMS_HitEnd_Prev = TMS_HitEnd;
              continue;
            }
            // thin layer is 3cm air + 1cm scint + 1.5cm steel = 5.5cm pitch
            // thick layer is 3cm air + 1cm scint + 4cm steel = 8cm pitch
            // if "gap" isn't a multiple of that, then something has gone off the rails and we want to stop
            double TMS_gap = TMS_HitStart.z()-TMS_HitStart_Prev.z();
            // correct for cosine of the incident angle
            double TMS_dist = (TMS_HitStart-TMS_HitStart_Prev).Mag();

            // If the gap between subsequent hits is less than the scintillator bar
            // We've had at least two hits in the same bar: skip these since the energy is already added
            if (TMS_gap < TMS_Thin_Layer_gap) {
              TMS_HitStart_Prev = TMS_HitStart;
              TMS_HitEnd_Prev = TMS_HitEnd;
              continue;
            }

            // If we're in the thin layer
            if (TMS_HitStart.z() < TMS_Thin_Layer_z_end) {
              int nplanes = TMS_gap/TMS_Thin_Layer_gap;
              if (nplanes > TMS_Max_Planes) {
                std::cerr << "Found a segment spanning more planes than the maximum!" << std::endl;
                std::cerr << "Diff: " << TMS_gap << " Gap: " << TMS_Thin_Layer_gap << " N_planes: " << nplanes << std::endl;

                std::cerr << "Previous hit start point: " << TMS_HitStart_Prev.Mag() << std::endl;
                std::cerr << "Current hit start point:  " << TMS_HitStart.Mag() << std::endl;
                std::cerr << "Event: " << ievt << std::endl;
                TMS_HitStart_Prev = TMS_HitStart;
                TMS_HitEnd_Prev = TMS_HitEnd;
                continue;
                //break;
              }
              TMS_total_planes += nplanes;
              // Have gone through 1cm scintillator, 1.5cm steel, 3cm air divided by cos(theta)
              TMS_TrackLength += (1.0*1.05+1.5*7.85)*nplanes*TMS_dist/TMS_gap;
              //std::cout << "Simply length of the track: " << (1+1.5+3)*nplanes*TMS_dist/TMS_gap << std::endl;
              //std::cout << "Full: " << (1.0*1.05+1.5*7.85)*nplanes*TMS_dist/TMS_gap << std::endl;
              // Without density
              //std::cout << "Scint only: " << (1.0*1.05)*nplanes*TMS_dist/TMS_gap << std::endl;

              // If we're in transition module, where the gap is 9.5
            } else if (TMS_HitStart.z()-TMS_Transition_Layer_z_end < 0.1) {
              //if (fabs(TMS_gap-TMS_Transition_Layer_gap) > 1.E-6) continue;
              if (fabs(TMS_gap-TMS_Transition_Layer_gap) > 1.E-6) break;
              int nplanes = TMS_gap/TMS_Transition_Layer_gap;
              // Increment
              TMS_total_planes += nplanes;
              // Have gone through 1cm scintillator, 4cm steel, 4.5cm air
              TMS_TrackLength += (1.0*1.05+4.0*7.85)*nplanes*TMS_dist/TMS_gap;
              //std::cout << "Simply length of the track: " << (7.85+1)*nplanes*TMS_dist/TMS_gap << std::endl;
              //std::cout << "Full: " << (1.0*1.05+1.5*7.85)*nplanes*TMS_dist/TMS_gap << std::endl;
              //std::cout << "Scint only: " << (1.0*1.05)*nplanes*TMS_dist/TMS_gap << std::endl;
              // In the thick region, where the gap is 8cm
            } else {
              int nplanes = TMS_gap/TMS_Thick_Layer_gap;
              //if (fabs(nplanes-TMS_gap/TMS_Thick_Layer_gap) > 1.E-6) continue;
              if (fabs(nplanes-TMS_gap/TMS_Thick_Layer_gap) > 1.E-6) break;
              if (nplanes > TMS_Max_Planes) {
                std::cerr << "Found a segment spanning more planes than the maximum!" << std::endl;
                std::cerr << "Diff: " << TMS_gap << " Gap: " << TMS_Thick_Layer_gap << " N_planes: " << nplanes << std::endl;

                std::cerr << "Previous hit start point: " << TMS_HitStart_Prev.Mag() << std::endl;
                std::cerr << "Current hit start point:  " << TMS_HitStart.Mag() << std::endl;
                std::cerr << "Event: " << ievt << std::endl;
                TMS_HitStart_Prev = TMS_HitStart;
                TMS_HitEnd_Prev = TMS_HitEnd;
                continue;
                //break;
              }
              TMS_total_planes += nplanes;
              // Have gone through 1cm scintillator, 4cm steel, 3cm air
              TMS_TrackLength += (1.0*1.05+4.0*7.85)*nplanes*TMS_dist/TMS_gap;
            }

            // Update the previous hit info
            TMS_HitStart_Prev = TMS_HitStart;
            TMS_HitEnd_Prev = TMS_HitEnd;
          }

          MuonDeath[0] = TMS_HitEnd_Prev.x();
          MuonDeath[1] = TMS_HitEnd_Prev.y();
          MuonDeath[2] = TMS_HitEnd_Prev.z();

          MuonScintLen = TMS_TrackLength;
          MuonScintEnergy = TMS_EnergyDeposit;
        }

      } // End loop over TG4HitSegmentDetectors (all the segments in various detectors in the event)
      // Add in the extra track length
      MuonScintLen += extra_trk_length;

      // Look at events which have both segments in the LAr and in the TMS
      //if (SegInLAr == false || SegInTMS == false) continue;

      tree_out->Fill();
    } // End loop over the vertices in the event
  } // End loop over all the events
  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << N_entries << " (" << Timer.RealTime()/N_entries << " s/event)" << std::endl;

  std::cout << "Writing to " << output->GetName() << "..." << std::endl;
  output->cd();
  tree_out->Write();
  output->Close();

  return true;
}

int main(int argc, char **argv) {
  if (argc != 2 && argc !=3) {
    std::cerr << "Need one argument: [EDepSim output file] [Output filename]" << std::endl;
    return -1;
  }

  std::string EDepSimFile = std::string(argv[1]);

  std::string OutputFile;
  // If only two arguments are given, 
  if (argc == 2) {
    std::string filename = std::string(argv[1]);
    OutputFile = filename.substr(0, filename.find(".root"));
    OutputFile += "_output.root";
  } else {
    OutputFile = std::string(argv[2]);
  }

  bool ok = dumpSSRITree(EDepSimFile, OutputFile);
  if (ok) return 0;
  else return -1;
}
