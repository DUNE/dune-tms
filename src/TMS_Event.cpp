#include "TMS_Event.h"

// Initialise the event counter to 0
int TMS_Event::EventCounter = 0;

TMS_Event::TMS_Event() {
  EventNumber = -999;
  SliceNumber = 0;
  SpillNumber = -999;
  nTrueTrajectories = -999;
  VertexIdOfMostEnergyInEvent = -999;
  LightWeight = true;
}

// Start the relatively tedious process of converting into TMS products!
// Can also use FillEvent = false to get a simple meta data extractor
TMS_Event::TMS_Event(TG4Event &event, bool FillEvent) {
  //std::cout<<"Making TMS event"<<std::endl;

  // Maybe make these class members
  // Keep false to process all events and all particles in events
  bool OnlyMuon = false;
  bool TMSOnly = false;
  bool TMSLArOnly = false;
  bool OnlyPrimary = true;
  bool LightWeight = TMS_Manager::GetInstance().Get_LightWeight_Truth();
  /*
  if (LightWeight) {
    OnlyMuon = true;
    TMSOnly = true;
    TMSLArOnly = true;
    OnlyPrimary = true;
  }
  */

  // Save down the event number
  EventNumber = EventCounter;
  SliceNumber = 0;
  SpillNumber = EventCounter;
  NSlices = 1; // By default there's at least one
  VertexIdOfMostEnergyInEvent = -999;

  // Check the integrity of the event
  //CheckIntegrity();

  int vtxcounter = 0;
  // Loop over the primary vertices
  for (TG4PrimaryVertexContainer::iterator it = event.Primaries.begin(); it != event.Primaries.end(); ++it) {
    //std::cout<<"for each event.Primaries "<<vtxcounter<<std::endl;

    TG4PrimaryVertex vtx = *it;
    Reaction = (*it).GetReaction();

    if (FillEvent) {
      // Primary particles in edep-sim are before any particle propagation happens
      // i.e. it's the particles out of the neutrino generation; don't save them
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;

      // Loop over the particles in the vertex and save them
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
        TMS_TrueParticle truepart = TMS_TrueParticle(particle, vtx);
        // Associate the particle with the position
        truepart.SetVertexID(vtxcounter);
        TMS_TruePrimaryParticles.emplace_back(truepart);
      }

      // Number of true trajectories
      nTrueTrajectories = event.Trajectories.size();
      // Now loop over the true trajectories (tracks) in the event
      for (TG4TrajectoryContainer::iterator jt = event.Trajectories.begin(); jt != event.Trajectories.end(); ++jt) {
        TG4Trajectory traj = *jt;

        // Only the muon if requested
        int PDGcode = traj.GetPDGCode();
        if (OnlyMuon && abs(PDGcode) != 13) continue;

        // Only from fundamental vertex if requested
        int ParentId = traj.GetParentId();
        if (OnlyPrimary && ParentId != -1) continue;

        // The id of this trajectory
        int TrackId = traj.GetTrackId();
        //std::cout << "PDG: " << PDGcode << " parentid: " << ParentId << " trackid: " << TrackId << " points: " << traj.Points.size() << std::endl;

        // Ignore particles that leave few hits, or gammas, if requested
        if (LightWeight && 
            //(traj.Points.size() < 3 || PDGcode == 22 || PDGcode == 2112)) continue;
            (PDGcode == 22 || PDGcode == 2112)) continue;
            //ParentId != -1) continue;

        // Is this the first time we encounter this particle in the trajectory point loop?
        bool firsttime = true;
        // Loop over the trajectory points of given true trajectory
        for (std::vector<TG4TrajectoryPoint>::iterator kt = traj.Points.begin(); kt != traj.Points.end(); kt++) {
          TG4TrajectoryPoint pt = *kt;

          // Check the point against the geometry
          TGeoNode *vol = TMS_Geom::GetInstance().GetGeometry()->FindNode(pt.GetPosition().X(), pt.GetPosition().Y(), pt.GetPosition().Z());

          // Very rarely but it does happen, the volume is null
          if (!vol) continue;
          std::string VolumeName = vol->GetName();

          // If asked to only look at LAr and TMS trajectories
          //
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

          // If firsttime is true and the above passes, this is the first time the particle enters any volume of interest, so create it
          if (firsttime) {
            // Can't set start momentum and position whe looping over the trajectory points, do this later
            //TMS_TrueParticle part(ParentId, TrackId, PDGcode, Momentum, Position);
            TMS_TrueParticle part(ParentId, TrackId, PDGcode);
            // Make the true particle that created this trajectory
            TMS_TrueParticles.push_back(std::move(part));
          } // End if (firsttime)

          // At this point we have a trajectory point that we are interested in, great!
          // Remember to fill this event with vertex information
          firsttime = false;

          // Now push back the position and momentum for the true particle at this trajectory point
          TLorentzVector Position = pt.GetPosition();
          TVector3 Momentum = pt.GetMomentum();

          // Might not want to save this truth information?
          // See G4ProcessType.hh, G4HaronicprocessType.hh, G4EmProcessSubType.hh
          int G4Process = pt.GetProcess();
          int G4Subprocess = pt.GetSubprocess();

          // Add the point
          TMS_TrueParticle *part = &(TMS_TrueParticles.back());
          part->AddPoint(Position, Momentum, G4Process, G4Subprocess);
        } // End loop over trajectory points

        // Save the birth and death points of trajectories that had a hit in a volume of interest
        if (!firsttime) {
          TG4TrajectoryPoint start = traj.Points.front();
          TG4TrajectoryPoint stop = traj.Points.back();

          // Get the TMS_TrueParticle corresponding to this particle
          TMS_TrueParticle *part = &(TMS_TrueParticles.back());

          TVector3 initialmom = start.GetMomentum();
          TLorentzVector initialpos = start.GetPosition();
          part->SetBirthMomentum(initialmom);
          part->SetBirthPosition(initialpos);

          TVector3 finalmom = stop.GetMomentum();
          TLorentzVector finalpos = stop.GetPosition();
          part->SetDeathMomentum(finalmom);
          part->SetDeathPosition(finalpos);
        } // End if (!firsttime) 
      } // End loop over the trajectories
      vtxcounter++;
    } // End if (FillEvent)
  } // End loop over the primary vertices, for (TG4PrimaryVertexContainer::iterator it
  
  // First create a mapping so we don't loop multiple times
  std::map<int, int> mapping_track_to_vertex_id;
  int vertex_index = 0;
  for (auto vertex : event.Primaries) {
    for (auto particle : vertex.Particles) {
      int track_id = particle.GetTrackId();
      mapping_track_to_vertex_id[track_id] = vertex_index;
    }
    vertex_index += 1;
  } 
  
  // Loop over each hit
  for (TG4HitSegmentDetectors::iterator jt = event.SegmentDetectors.begin(); jt != event.SegmentDetectors.end(); ++jt) {
    // Only look at TMS hits
    std::string DetString = (*jt).first;

    // Skip hits outside of the TMS if running lightweight
    if (TMSOnly && DetString != TMS_Const::TMS_EDepSim_VolumeName) continue;

    TG4HitSegmentContainer tms_hits = (*jt).second;
    for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
      TG4HitSegment edep_hit = *kt;
      int track_id = edep_hit.GetPrimaryId();
      int vertex_id = mapping_track_to_vertex_id[track_id];
      TMS_Hit hit = TMS_Hit(edep_hit, vertex_id);
      TMS_Hits.push_back(std::move(hit));
      
      // todo, maybe skip for michel electrons or late neutrons
      TrueVisibleEnergyPerVertex[hit.GetTrueHit().GetVertexId()] += hit.GetTrueHit().GetE();

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
    } // End for (TG4HitSegmentContainer::iterator kt
  } // End loop over each hit, for (TG4HitSegmentDetectors::iterator jt
  
  // Merge hits that happened in the same scintillator strip and within the same readout time window
  // This is a simulation cleanup step, not reconstruction
  MergeCoincidentHits();
  // Now apply optical and timing models
  ApplyReconstructionEffects();

  EventCounter++;
}

TMS_Event::TMS_Event(TMS_Event &event, int slice) {
  // Create an event from a slice of another event
  TMS_Hits = event.GetHits(slice);
  SliceNumber = slice;
  SpillNumber = event.SpillNumber;
  // Todo, did I copy everything
  // Update event counter if slice != 0, and keep old event number for slice 0.
  if (slice != 0) {
    EventNumber = EventCounter;
    EventCounter++;
  }
  else {
    EventNumber = event.EventNumber;
  }
}

void TMS_Event::MergeCoincidentHits() {
  std::sort(TMS_Hits.begin(), TMS_Hits.end(), TMS_Hit::SortByZThenT);
  
  // Loop over the original hits
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); 
      it != TMS_Hits.end(); it++) {
    if ((*it).GetPedSup()) continue; // Skip hits which are already removed
    // Maybe this hit has already been counted
    double z = (*it).GetZ();
    double y = (*it).GetNotZ();
    //double e = hit.GetE();
    double t = (*it).GetT();
    
    //Strip out hits that are outside the actual volume 
    // This is probably some bug in the geometry that sometimes gives hits in the z=30k mm (i.e. 10m downstream of the end of the TMS)
    // TODO figure out why these happen
    if (z > TMS_Const::TMS_End[2] || z < TMS_Const::TMS_Start[2]) (*it).SetPedSup(true);

    // Look ahead to find duplicates, but stop when z != z2
    std::vector<std::vector<TMS_Hit>::iterator> duplicates;
    for (std::vector<TMS_Hit>::iterator jt = it+1; jt != TMS_Hits.end() && jt->GetZ() == z; ++jt) {

      TMS_Hit hit2 = *(jt);
      if (hit2.GetPedSup()) continue; // Skip hits which are already removed
      double z2 = hit2.GetZ();
      double y2 = hit2.GetNotZ();
      double e2 = hit2.GetE();
      double t2 = hit2.GetT();

      // Merge
      if (z == z2 && y == y2 && fabs(t2-t) < TMS_Const::TMS_TimeThreshold) {
        (*it).SetE((*it).GetE()+e2);
        (*it).SetT(std::min(t, t2));
        duplicates.push_back(jt);
      }
    }
    // Now flag the duplicates for removal
    for (int i = duplicates.size() - 1; i >= 0; i--) {
      auto hit_to_erase = duplicates[i];
      hit_to_erase->SetPedSup(true);
    }
  }
  // Now erase all hits that are set as ped supped
  std::vector<TMS_Hit> remaining_hits;
  std::vector<TMS_Hit> deleted_hits;
  for (auto hit : TMS_Hits) {
    if (!hit.GetPedSup()) remaining_hits.push_back(hit);
    else deleted_hits.push_back(hit);
  }
  TMS_Hits.clear();
  for (auto hit : remaining_hits) TMS_Hits.push_back(hit);
  deleted_hits.erase(deleted_hits.begin(), deleted_hits.end());
}

void TMS_Event::SimulateOpticalModel() {
  
}

void TMS_Event::SimulatePedestalSubtraction() {
  // Don't actually remove the hits because they may be relevant for other processes
  // Loop over the hits
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); it++) {
     auto hit = (*it);
     double hit_energy = hit.GetE();
     if (hit_energy < TMS_Const::TMS_EnThres) {
       hit.SetPedSup(true);
     }   
  }
}

void TMS_Event::SimulateTimingModel() {
  // List of timing effects to simulate:
  // Random electronic timing noise
  // Deadtime
  // Time skew from first PE to hit sensor
  // Optical fiber length delays (corrected to strip center)
  // Timing effects from random noise, cross talk, afterpulsing
}

void TMS_Event::ApplyReconstructionEffects() {
  // Simulate an optical model 
  SimulateOpticalModel();
  // Simulate pedestal subtraction where any hit under TMS_Const::TMS_EnThres MeV is removed
  SimulatePedestalSubtraction();
  // Simulate a timing model
  SimulateTimingModel();
}

const std::vector<TMS_Hit> TMS_Event::GetHits(int slice, bool include_ped_sup) {
  std::vector<TMS_Hit> out;
  for (auto hit : TMS_Hits) {
    if (!hit.GetPedSup() || include_ped_sup) {
      int slice_number = hit.GetSlice();
      //std::cout<<"Slice number for hit: "<<slice_number<<std::endl;
      if (slice_number == slice || slice < 0) out.push_back(hit);
    }
  }
  return out;
}

// Add a separate event to this event
// Handy for making hacked overlays
void TMS_Event::AddEvent(TMS_Event &Other_Event) {
  std::cout << "Adding event " << Other_Event.GetEventNumber() << " to event " << GetEventNumber() << std::endl;

  // Get the other hits
  std::vector<TMS_Hit> other_hits = Other_Event.GetHits(-1, true);

  // And use them to expand on the original hits in the event
  for (auto &hit: other_hits) {
    TMS_Hits.emplace_back(std::move(hit));
  }

  // Do the same for the true particles
  std::vector<TMS_TrueParticle> other_truepart = Other_Event.GetTrueParticles();
  for (auto &part: other_truepart) {
    TMS_TrueParticles.emplace_back(std::move(part));
  }

}

// For now just fill the true neutrino
// But shows how you can easily make a vector of rootracker particles for the TMS_Event to carry along
void TMS_Event::FillTruthFromGRooTracker(int pdg[__EDEP_SIM_MAX_PART__], double p4[__EDEP_SIM_MAX_PART__][4]) {
  TrueNeutrino.first.SetX(p4[0][0]);
  TrueNeutrino.first.SetY(p4[0][1]);
  TrueNeutrino.first.SetZ(p4[0][2]);
  TrueNeutrino.first.SetT(p4[0][3]);
  TrueNeutrino.second = pdg[0];
  
}

void TMS_Event::FillAdditionalTruthFromGRooTracker(double x4[__EDEP_SIM_MAX_PART__][4]) {

  TrueNeutrinoPosition.SetX(x4[0][0]);
  TrueNeutrinoPosition.SetY(x4[0][1]);
  TrueNeutrinoPosition.SetZ(x4[0][2]);
  TrueNeutrinoPosition.SetT(x4[0][3]);
}

void TMS_Event::FillTrueLeptonInfo(int pdg, TLorentzVector position, TLorentzVector momentum) {
  TrueLeptonPDG = pdg;
  TrueLeptonPosition = position;
  TrueLeptonMomentum = momentum;
}

int TMS_Event::GetVertexIdOfMostVisibleEnergy() {
  // Return early if we've already calculated it
  if (VertexIdOfMostEnergyInEvent >= 0) return VertexIdOfMostEnergyInEvent;

  // Reset the map
  TrueVisibleEnergyPerVertex.clear();
  // First find how much energy is in each variable
  for (auto hit : TMS_Hits) {
    int vertex_id = hit.GetTrueHit().GetVertexId();
    // todo, true or reco energy?
    double energy = hit.GetTrueHit().GetE();
    TrueVisibleEnergyPerVertex[vertex_id] += energy;
  }
  
  // Now find the most energetic vertex
  double max = -1e9;
  int max_vertex_id = -999;
  double total_energy = 0;
  for (auto it : TrueVisibleEnergyPerVertex) {
    double vertex = it.first;
    double energy = it.second;
    //std::cout<<"Vertex "<<vertex<<" has E: "<<energy<<std::endl;
    if (energy > max) { max = energy; max_vertex_id = vertex; }
    total_energy += energy;
  }
  VertexIdOfMostEnergyInEvent = max_vertex_id;
  VisibleEnergyFromVertexInSlice = max;
  VisibleEnergyFromOtherVerticesInSlice = total_energy - max;
  
  return VertexIdOfMostEnergyInEvent;
}

void TMS_Event::FillTruthInformation(TG4Event &event) {
  // First create a mapping so we don't loop twice
  /*std::map<int, int> mapping;
  int vertex_index = 0;
  for (auto vertex : event.Primaries) {
    for (auto particle : vertex.Particles) {
      int track_id = particle.GetTrackId();
      mapping[track_id] = vertex_index;
    }
    vertex_index += 1;
  } 
  
  for (size_t i = 0; i < TMS_Hits.size(); i++) {
    TMS_Hit hit = TMS_Hits[i];
    int primary_id = hit.GetTrueHit().GetPrimaryId();
    auto traj = event.Trajectories[primary_id];
    int parent_id = traj.GetParentId();
    // Keep going up the tree until you find the primary particle
    while (parent_id >= 0) {
      traj = event.Trajectories[parent_id];
      primary_id = parent_id;
      parent_id = traj.GetParentId();
    }
    int vertex_id = mapping[primary_id];
    std::cout<<"primary_id: "<<primary_id<<", vertex id: "<<vertex_id<<std::endl;
    if (vertex_id >= 0) {
      hit.GetTrueHit().SetVertexId(vertex_id);
    }
  }*/
  // TODO erase function
  event.Primaries.begin();
}

std::pair<double, double> TMS_Event::GetEventTimeRange() {
  double min_time = 1e9;
  double max_time = -1e9;
  for (auto hit : TMS_Hits) {
    min_time = std::min(min_time, hit.GetT());
    max_time = std::max(max_time, hit.GetT());
  }
  return std::make_pair(min_time, max_time);
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
  std::cout << "  Vertex ID of most energy: " << VertexIdOfMostEnergyInEvent << std::endl;
  std::cout << "  Visible energy in slice: " << VisibleEnergyFromVertexInSlice << std::endl;
  std::cout << "  Total visible energy: " << TotalVisibleEnergyFromVertex << std::endl;
  std::cout << "  Other visible energy: " << VisibleEnergyFromOtherVerticesInSlice << std::endl;

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

double TMS_Event::GetMuonTrueKE() {
  std::vector<TMS_TrueParticle> TrueParticles = GetTrueParticles();
  double HighestKE = -999.99;
  for (auto it = TrueParticles.begin(); it != TrueParticles.end(); ++it) {
    // Only save muon info for now
    if (abs((*it).GetPDG()) != 13) continue;
    // Also make sure it's a fundamental muon
    if ((*it).GetParent() != -1) continue;
    TVector3 mom = (*it).GetBirthMomentum();
    double E = (*it).GetBirthEnergy();
    // Get KE (E - m)
    double mass = sqrt(E*E-mom.Mag2());
    double KE = E-mass;
    if (KE > HighestKE) HighestKE = KE;
  }
  return HighestKE;
}

double TMS_Event::GetMuonTrueTrackLength() {
  std::vector<TMS_TrueParticle> TrueParticles = GetTrueParticles();
  double total = 0;
  for (auto it = TrueParticles.begin(); it != TrueParticles.end(); ++it) {
    // Only save muon info for now
    if (abs((*it).GetPDG()) != 13) continue;
    // Also make sure it's a fundamental muon
    if ((*it).GetParent() != -1) continue;

    std::vector<TLorentzVector> pos = (*it).GetPositionPoints();
    int num = 0;
    for (auto pnt = pos.begin(); pnt != pos.end(); ++pnt, ++num) {
      auto nextpnt = *(pnt+1);
      TVector3 point1((*pnt).X(), -200, (*pnt).Z());
      TVector3 point2(nextpnt.X(), -200, nextpnt.Z());
      //point1.Print();
      //point2.Print();
      if ((point2-point1).Mag() > 100) {
        //std::cout << "moving on" << std::endl;
        continue;
      }
      double tracklength = TMS_Geom::GetInstance().GetTrackLength(point1, point2);
      total += tracklength;
      //std::cout << "point " << num << std::endl;
      //std::cout << "total: " << total << std::endl;
      //std::cout << "tracklength: " << tracklength << std::endl;
    }
  }
  return total;
}
