#include "TMS_Event.h"

// Initialise the event counter to 0
int TMS_Event::EventCounter = 0;

// Start the relatively tedious process of converting into TMS products!
// Can also use FillEvent = false to get a simple meta data extractor
TMS_Event::TMS_Event(TG4Event &event, bool FillEvent) {

  // Maybe make these class members
  bool OnlyMuon = false;
  bool TMSOnly = false;
  bool TMSLArOnly = false;
  bool OnlyPrimary = false;
  bool LightWeight = TMS_Manager::GetInstance().Get_LightWeight_Truth();

  // Save down the event number
  EventNumber = EventCounter;

  // Check the integrity of the event
  //CheckIntegrity();

  int vtxcounter = 0;
  // Loop over the primary vertices
  for (TG4PrimaryVertexContainer::iterator it = event.Primaries.begin(); it != event.Primaries.end(); ++it) {

    TG4PrimaryVertex vtx = *it;
    Reaction = (*it).GetReaction();

    if (FillEvent) {
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;

      // Primary particles in edep-sim are before any particle propagation happens
      // i.e. it's the particles out of the neutrino generation; don't save them

      // Loop over the particles in the vertex and save them
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
        TMS_TrueParticle truepart = TMS_TrueParticle(particle, vtx);
        // Associate the particle with the position
        truepart.SetVertexID(vtxcounter);
        TMS_TruePrimaryParticles.emplace_back(truepart);
      }

      nTrueTrajectories = event.Trajectories.size();

      // Now loop over the trajectories (tracks) in the event
      for (TG4TrajectoryContainer::iterator jt = event.Trajectories.begin(); jt != event.Trajectories.end(); ++jt) {
        TG4Trajectory traj = *jt;

        // Only the muon
        int PDGcode = traj.GetPDGCode();
        if (OnlyMuon && abs(PDGcode) != 13) continue;

        // Only from fundamental vertex
        int ParentId = traj.GetParentId();
        if (OnlyPrimary && ParentId != -1) continue;

        // The id of this trajectory
        int TrackId = traj.GetTrackId();
        //std::cout << "PDG: " << PDGcode << " parentid: " << ParentId << " trackid: " << TrackId << " points: " << traj.Points.size() << std::endl;

        // Ignore particles that leave few hits, or gammas
        if (LightWeight && 
            (traj.Points.size() < 3 || PDGcode == 22 || PDGcode == 2112)) continue;
            //ParentId != -1) continue;

        // Skip this trajectory if it has no parts of it in TMS or LAr
        bool skip = true;

        // Save down the trajectory points of the true particles
        for (std::vector<TG4TrajectoryPoint>::iterator kt = traj.Points.begin(); kt != traj.Points.end(); kt++) {
          TG4TrajectoryPoint pt = *kt;

          // Check the point against the geometry
          TGeoNode *vol = TMS_Geom::GetInstance().GetGeometry()->FindNode(pt.GetPosition().X(), pt.GetPosition().Y(), pt.GetPosition().Z());

          // Very rarely but it does happen, the volume is null
          if (!vol) continue;
          std::string VolumeName = vol->GetName();

          // Only look at LAr and TMS hits
          if (TMSOnly || TMSLArOnly) {
            // Check the TMS volume first
            if (VolumeName.find(TMS_Const::TMS_VolumeName) == std::string::npos && 
                VolumeName.find(TMS_Const::TMS_ModuleLayerName) == std::string::npos &&
                VolumeName.find(TMS_Const::TMS_EDepSim_VolumeName) == std::string::npos) continue;

            // check the LAr volume
            if (TMSLArOnly) {
              if (VolumeName.find(TMS_Const::LAr_ActiveName) == std::string::npos) continue;
            }
          }

          // If skip is true and the above passes, this is the first time the particle enters volume of interest, so create it
          if (skip) {
            //TMS_TrueParticle part(ParentId, TrackId, PDGcode, Momentum, Position);
            TMS_TrueParticle part(ParentId, TrackId, PDGcode);
            // Make the true particle that created this trajectory
            TMS_TrueParticles.push_back(std::move(part));
          }

          // At this point we have a trajectory point in the TMS, great!
          // Remember to fill this event with vertex information
          skip = false;

          // Now push back the position and momentum for the true particle at this trajectory point
          TLorentzVector Position = pt.GetPosition();
          TVector3 Momentum = pt.GetMomentum();
          // Add the point
          TMS_TrueParticle *part = &(TMS_TrueParticles.back());
          part->AddPoint(Position, Momentum);
        } // End loop over trajectory point

        // Save the birth and death points of trajectory that had a hit in TMS or LAr
        if (!skip) {
          TG4TrajectoryPoint start = traj.Points.front();
          TG4TrajectoryPoint stop = traj.Points.back();

          TMS_TrueParticle *part = &(TMS_TrueParticles.back());

          TVector3 initialmom = start.GetMomentum();
          TLorentzVector initialpos = start.GetPosition();
          part->SetBirthMomentum(initialmom);
          part->SetBirthPosition(initialpos);

          TVector3 finalmom = stop.GetMomentum();
          TLorentzVector finalpos = stop.GetPosition();
          part->SetDeathMomentum(finalmom);
          part->SetDeathPosition(finalpos);
        }
      }

      // Loop over each hit
      for (TG4HitSegmentDetectors::iterator jt = event.SegmentDetectors.begin(); jt != event.SegmentDetectors.end(); ++jt) {
        // Only look at TMS hits
        std::string DetString = (*jt).first;
        if (DetString != TMS_Const::TMS_EDepSim_VolumeName) continue;

        TG4HitSegmentContainer tms_hits = (*jt).second;
        for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
          TG4HitSegment edep_hit = *kt;
          TMS_Hit hit = TMS_Hit(edep_hit);
          TMS_Hits.push_back(std::move(hit));


          // Loop through the True Particle list and associate
          /*
          // Now associate the hits with the muon
          int PrimaryId = edep_hit.GetPrimaryId();
          for (auto &TrueParticle : TMS_TrueParticles) {
          // Check the primary ID
          if (TrueParticle.GetTrackId() != PrimaryId) continue;
          TLorentzVector Position = (edep_hit.GetStop()+edep_hit.GetStart());
          Position *= 0.5;
          TrueParticle.AddPoint(Position);
          }
          */
        }
      }
      vtxcounter++;
    }
  }

  EventCounter++;
}

void TMS_Event::FillTruthFromGRooTracker(int pdg[__EDEP_SIM_MAX_PART__], double p4[__EDEP_SIM_MAX_PART__][4]) {
  TrueNeutrino.first.SetX(p4[0][0]);
  TrueNeutrino.first.SetY(p4[0][1]);
  TrueNeutrino.first.SetZ(p4[0][2]);
  TrueNeutrino.first.SetT(p4[0][3]);
  TrueNeutrino.second = pdg[0];
}

void TMS_Event::Print() {
  std::cout << std::endl;
  std::cout << "*** " << std::endl;
  std::cout << "Printing TMS_Event class from "  << __FILE__ << std::endl;
  std::cout << "  Using geometry: " << TMS_Geom::GetInstance().GetGeometry()->GetName() << ", " << TMS_Geom::GetInstance().GetGeometry()->GetTitle() << std::endl;
  std::cout << "  From: " << TMS_Geom::GetInstance().GetFileName() << std::endl;
  std::cout << "  Initial state neutrino from gRooTracker: " << std::endl;
  std::cout << "  PDG: " << TrueNeutrino.second << " (px, py, pz, E) = (" << TrueNeutrino.first.X() << ", " << TrueNeutrino.first.Y() << ", " << TrueNeutrino.first.Z() << ", " << TrueNeutrino.first.T() << ")" << std::endl;
  std::cout << "  N True primary particles: " << TMS_TruePrimaryParticles.size() << std::endl;
  std::cout << "  N True filtered trajectories: " << TMS_TrueParticles.size() << std::endl;
  std::cout << "  N True unfiltered trajectories: " << nTrueTrajectories << std::endl;
  std::cout << "  N Hits: " << TMS_Hits.size() << std::endl;
  std::cout << "  IsEmpty: " << IsEmpty() << std::endl;

  std::cout << "Printing primary particle stack: " << std::endl;
  int PartCount = 0;
  for (auto it = TMS_TruePrimaryParticles.begin(); it != TMS_TruePrimaryParticles.end(); ++it, ++PartCount) {
    std::cout << "Particle " << PartCount << std::endl;
    (*it).Print();
  }

  PartCount = 0;
  std::cout << "Printing trajectory particle stack: " << std::endl;
  for (std::vector<TMS_TrueParticle>::iterator it = TMS_TrueParticles.begin(); it != TMS_TrueParticles.end(); ++it, ++PartCount) {
    std::cout << "Particle " << PartCount << std::endl;
    (*it).Print();
  }

  int HitCount = 0;
  std::cout << "Printing hit stack: " << std::endl;
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it, ++HitCount) {
    std::cout << "Hit "  << HitCount << std::endl;
    (*it).Print();
  }
}

