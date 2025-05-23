#include "TMS_Event.h"
#include "TMS_Readout_Manager.h"
#include "TDatabasePDG.h"
#include <random>

// Initialise the event counter to 0
int TMS_Event::EventCounter = 0;

TMS_Event::TMS_Event() {
  EventNumber = -999;
  SliceNumber = 0;
  SpillNumber = -999;
  nTrueTrajectories = -999;
  nVertices = -999;
  VertexIdOfMostEnergyInEvent = -999;
  LightWeight = true;
}

static bool TMS_TrueParticle_NotWorthSaving(TMS_TrueParticle tp) {
  if (tp.GetTrueVisibleEnergy() == 0 && !tp.IsPrimary()) return true;
  // Don't worry about really low visible energy
  if (tp.GetTrueVisibleEnergy() < 0.5 && !tp.IsPrimary()) return true;
  else return false;
};

void TMS_Event::ProcessTG4Event(TG4Event &event, bool FillEvent) {

  TDatabasePDG *database = TDatabasePDG::Instance();
  
  // Maybe make these class members
  // Keep false to process all events and all particles in events
  bool OnlyMuon = false;
  bool TMSOnly = false;
  bool TMSLArOnly = false;
  bool OnlyPrimary = false;
  bool OnlyPrimaryOrInteresting = false;
  bool LightWeight = TMS_Manager::GetInstance().Get_LightWeight_Truth();
  
  int nPrimary = 0;
  int nInteresting = 0;
  int nTotal = 0;
  int nCharged = 0;
  int nHighMomentum = 0;
  int nChargedAndLowMomentum = 0;
  RunNumber = event.RunId;
  int current_vertexid = event.EventId;
  // Nersc jobs have 1 primary vertex per entry, whereas fermigrid jobs have many, but don't use the spill builder.
  // So they're not affected by the https://github.com/DUNE/2x2_sim/issues/54 bug
  // todo: Make both GetInteractionNumber when bug is fixed
  bool use_GetInteractionNumber = false;
  if (event.Primaries.size() > 1) use_GetInteractionNumber = true; 
  // Loop over the primary vertices
  for (TG4PrimaryVertexContainer::iterator it = event.Primaries.begin(); it != event.Primaries.end(); ++it) {

    TG4PrimaryVertex vtx = *it;
    Reaction = (*it).GetReaction();
    
    // Interaction number is off-by-one in recent microprod files, so set it manually
    // See https://github.com/DUNE/2x2_sim/issues/61
    vtx.InteractionNumber = current_vertexid;
    // Ideally we'd do it like this, but it's not supported by the spill builder
    // See https://github.com/DUNE/2x2_sim/issues/54
    if (use_GetInteractionNumber)
      current_vertexid = vtx.GetInteractionNumber();
    if (current_vertexid < 0) {
      std::cout<<"Fatal: Got a current_vertexid < 0 in TMS_Event: "<<current_vertexid<<std::endl;
      throw std::runtime_error("Fatal: Get a vertex id < 0");
    }
    Reactions[current_vertexid] = Reaction;
      
    Vtx_Info vtx_info;
    vtx_info.reaction = Reaction;
    vtx_info.vtx_id = current_vertexid;
    // Had issues with lorentz vectors before so best make a copy
    vtx_info.SetVtx(TLorentzVector(vtx.GetPosition().X(), vtx.GetPosition().Y(), vtx.GetPosition().Z(), vtx.GetPosition().T()));
    
    info_about_vtx[current_vertexid] = vtx_info;

    if (FillEvent) {
      // Primary particles in edep-sim are before any particle propagation happens
      // i.e. it's the particles out of the neutrino generation; don't save them
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;

      // Loop over the particles in the vertex and save them
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
        TMS_TrueParticle truepart = TMS_TrueParticle(particle, vtx);
        TMS_TruePrimaryParticles.emplace_back(truepart);
        
        if (current_vertexid != truepart.GetVertexID()) {
          std::cout<<"Fatal: TMS_TrueParticle's vertex id was set incorrectly in true primary particle list" \
                     " and doesn't match the current id: true part vtx id: ";
          std::cout<<truepart.GetVertexID()<<" vs current: "<<current_vertexid<<std::endl;
          throw std::runtime_error("Fatal: TMS_TrueParticle's vertex id was set incorrectly");
        }
      }

      // Number of true trajectories
      nTrueTrajectories = event.Trajectories.size();
      // Now loop over the true trajectories (tracks) in the event
      for (TG4TrajectoryContainer::iterator jt = event.Trajectories.begin(); jt != event.Trajectories.end(); ++jt) {
        nTotal += 1;
        TG4Trajectory traj = *jt;

        // Only the muon if requested
        int PDGcode = traj.GetPDGCode();
        if (OnlyMuon && abs(PDGcode) != 13) continue;

        // Only from fundamental vertex if requested
        int ParentId = traj.GetParentId();
        if (OnlyPrimary && ParentId != -1) continue;
        bool isPrimary = ParentId == -1;

        // The id of this trajectory
        int TrackId = traj.GetTrackId();
        //std::cout << "PDG: " << PDGcode << " parentid: " << ParentId << " trackid: " << TrackId << " points: " << traj.Points.size() << std::endl;

        // Ignore particles that leave few hits, or gammas, if requested
        if (LightWeight && 
            //(traj.Points.size() < 3 || PDGcode == 22 || PDGcode == 2112)) continue;
            (PDGcode == 22 || PDGcode == 2112)) continue;
            //ParentId != -1) continue;

        bool isCharged = false;
        bool isHighMomentum = false;
        if (PDGcode > 1000000000) {
          // Numbers above 1000000000 are nuclei, and so aren't in the database
          // They are charged though, but unlikely to have enough momentum to go far
          isCharged = true;
        } else {
          auto particle = database->GetParticle(PDGcode);
          if (!particle) {
            std::cout<<"Warning: Couldn't find pdg code "<<PDGcode<<" in pdg database"<<std::endl;
          }
          else {
            // Check if it's neutral or not
            isCharged = std::abs(particle->Charge()) > 0.2;
          }
        }
        if (traj.Points.size() > 0) {
          TVector3 initial_momentum = traj.Points[0].GetMomentum();
          if (initial_momentum.Mag() > 5) isHighMomentum = true;
          //if (isCharged && isHighMomentum && !isPrimary) {
          //    std::cout<<"Found interesting non-primary particle "<<PDGcode<<", momentum="<<initial_momentum.Mag()<<std::endl;
          //}
        }
        bool isInteresting = isHighMomentum && isCharged; 

        if (isPrimary) nPrimary += 1;
        if (!isPrimary && isInteresting) nInteresting += 1;
        if (!isPrimary) {
          if (isCharged) nCharged += 1;
          if (isHighMomentum) nHighMomentum += 1;
          if (isCharged && !isHighMomentum) nChargedAndLowMomentum += 1;
        }

        // Skip if not interesting and not primary
        if (OnlyPrimaryOrInteresting) {
          if ((!isPrimary) && (!isInteresting)) continue;
        }

        // Is this the first time we encounter this particle in the trajectory point loop?
        bool firsttime = true;
        // Loop over the trajectory points of given true trajectory
        for (std::vector<TG4TrajectoryPoint>::iterator kt = traj.Points.begin(); kt != traj.Points.end(); kt++) {
          TG4TrajectoryPoint pt = *kt;

          // Check the point against the geometry
          TGeoNode *vol = TMS_Geom::GetInstance().FindNode(pt.GetPosition().X(), pt.GetPosition().Y(), pt.GetPosition().Z());

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
            TMS_TrueParticle part(ParentId, TrackId, PDGcode, current_vertexid);
            // Make the true particle that created this trajectory
            TMS_TrueParticles.push_back(std::move(part));
        
            if (current_vertexid != part.GetVertexID()) {
              std::cout<<"Fatal: TMS_TrueParticle's vertex id was set incorrectly in all particle list " \
                         "and doesn't match the current id: true part vtx id: ";
              std::cout<<part.GetVertexID()<<" vs current: "<<current_vertexid<<std::endl;
              throw std::runtime_error("Fatal: TMS_TrueParticle's vertex id was set incorrectly in all particle list");
            }
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
      nVertices++;
    } // End if (FillEvent)
  } // End loop over the primary vertices, for (TG4PrimaryVertexContainer::iterator it
  
  
  //std::cout<<"N total: "<<nTotal<<", N Primary: "<<nPrimary<<", N Interesting: "<<nInteresting<<", N charged: "<<nCharged<<", N high P: "<<nHighMomentum<<", N charged and low P: "<<nChargedAndLowMomentum<<", n TMS_TruePrimaryParticles: "<<TMS_TruePrimaryParticles.size()<<std::endl;

  // First create a mapping so we don't loop multiple times
  std::map<int, int> mapping_track_to_vertex_id;
  int vertex_index = event.EventId;
  for (auto vertex : event.Primaries) {
    // Ideally we'd do it like this for nersc files, but it's not supported by the spill builder
    // See https://github.com/DUNE/2x2_sim/issues/54
    if (use_GetInteractionNumber)
      vertex_index = vertex.GetInteractionNumber();
    for (auto particle : vertex.Particles) {
      int track_id = particle.GetTrackId();
      mapping_track_to_vertex_id[track_id] = vertex_index;
    }
    for (auto traj : event.Trajectories) {
      int track_id = traj.GetTrackId();
      mapping_track_to_vertex_id[track_id] = vertex_index;
    }
  }

  std::map<int, TMS_TrueParticle*> mapping_track_to_true_particle;
  for (auto& tp : TMS_TrueParticles) {
    int key = tp.GetVertexID() * 100000 + tp.GetTrackId();
    mapping_track_to_true_particle[key] = &tp;
  }
  
  std::map<std::tuple<int, int, int, int>, size_t> map_pos_nontms_hits;

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
      int vertex_id = -999;
      auto value = mapping_track_to_vertex_id.find(track_id);
      if (value == mapping_track_to_vertex_id.end()) {
        std::cout<<"WARNING: Didn't find track id in mapping_track_to_vertex_id! track_id = "<<track_id<<", mapping_track_to_vertex_id.size() = "<<mapping_track_to_vertex_id.size()<<", this shouldn't happen anymore\n\n\n"<<std::endl;
      }
      else vertex_id = value->second;
      TMS_Hit hit = TMS_Hit(edep_hit, vertex_id);
      int barnum = hit.GetBarNumber();
      // Only add if within the TMS
      // Can't use x,y or z because geometry might change. But we know things aren't set if there's no bar number
      if (barnum >= 0) {
        auto t = hit.GetAdjustableTrueHit();
        for (size_t i = 0; i < t.GetNTrueParticles(); i++) {
          int key = t.GetVertexIds(i) * 100000 + t.GetPrimaryIds(i);
          if (mapping_track_to_true_particle.find(key) != mapping_track_to_true_particle.end()) {
            // Now set info
            auto tp = mapping_track_to_true_particle[key];
            if (tp->IsLeptonic()) t.SetEnergyLeptonic(i);
          }
        }
        SaveKeyVertexInfo(t);
        TMS_Hits.push_back(std::move(hit));

        // todo, maybe skip for michel electrons or late neutrons
        for (size_t i = 0; i < hit.GetTrueHit().GetNTrueParticles(); i++) {
          TrueVisibleEnergyPerVertex[hit.GetTrueHit().GetVertexIds(i)] += hit.GetTrueHit().GetEnergyShare(i);
          TrueVisibleEnergyPerParticle[hit.GetTrueHit().GetVertexIds(i) * 100000 + hit.GetTrueHit().GetPrimaryIds(i)] += hit.GetTrueHit().GetEnergyShare(i);
        }
      }
      else if (DetString.find(TMS_Const::LAr_ActiveName) != std::string::npos) {
        // Only care about LAr active volume
        // We only need it for truth info so just save truth info
        TMS_TrueHit t(edep_hit, vertex_id);
        for (size_t i = 0; i < t.GetNTrueParticles(); i++) {
          int key = t.GetVertexIds(i) * 100000 + t.GetPrimaryIds(i);
          auto itp = mapping_track_to_true_particle.find(key);
          if (itp != mapping_track_to_true_particle.end()) {
            // Now set info
            auto tp = itp->second;
            if (tp->IsLeptonic()) t.SetEnergyLeptonic(i);
          }
        }
        
        SaveKeyVertexInfo(t);
        
        double divide = 10.0;
        auto poskey = std::tuple((int) (t.GetX() / divide), (int) (t.GetY() / divide), (int) (t.GetZ() / divide), t.GetVertexIds(0));
        if (map_pos_nontms_hits.find(poskey) != map_pos_nontms_hits.end()) {
          // Already exists, merge with existing
          auto& merge_with_me = NonTMS_Hits[map_pos_nontms_hits[poskey]];
          merge_with_me.MergeWith(t);
        }
        else {
          // Doesn't exist, add to list and map
          NonTMS_Hits.push_back(t);
          map_pos_nontms_hits[poskey] = NonTMS_Hits.size() - 1;
        }
      }
    } // End for (TG4HitSegmentContainer::iterator kt
  } // End loop over each hit, for (TG4HitSegmentDetectors::iterator jt
  bool OnlyPrimaryOrVisibleEnergy = true;
  
  // Now update truth info per particle
  for (size_t i = 0; i < TMS_TrueParticles.size(); i++) {
    double energy = 0;
    // If it's not in the map, don't create it
    int key = TMS_TrueParticles[i].GetVertexID() * 100000 + TMS_TrueParticles[i].GetTrackId();
    auto it = TrueVisibleEnergyPerParticle.find(key);
    if (it != TrueVisibleEnergyPerParticle.end()) {
      energy = it->second;
    }
    TMS_TrueParticles[i].SetTrueVisibleEnergy(energy, false);
  }
  nTrueForgottenParticles = -1;
  if (OnlyPrimaryOrVisibleEnergy) {
    size_t initial = TMS_TrueParticles.size();
    TMS_TrueParticles.erase(std::remove_if(TMS_TrueParticles.begin(), 
                            TMS_TrueParticles.end(), 
                            TMS_TrueParticle_NotWorthSaving), 
                            TMS_TrueParticles.end());
    size_t end = TMS_TrueParticles.size();
    nTrueForgottenParticles = initial - end;
  }
}

// Start the relatively tedious process of converting into TMS products!
// Can also use FillEvent = false to get a simple meta data extractor
TMS_Event::TMS_Event(TG4Event event, bool FillEvent) {
  //std::cout<<"Making TMS event"<<std::endl;

  // Save down the event number
  EventNumber = EventCounter;
  generator = std::default_random_engine(7890 + EventNumber); 
  SliceNumber = 0;
  SpillNumber = EventCounter;
  NSlices = 1; // By default there's at least one
  VertexIdOfMostEnergyInEvent = -999;
  nVertices = 0;

  // Check the integrity of the event
  //CheckIntegrity();

  ProcessTG4Event(event, FillEvent);
  
  // Now apply optical and timing models
  //ApplyReconstructionEffects();
  // TODO figure out why SimulateOpticalModel/MergeCoincidentHits are needed here to avoid crash
  // Simulate an optical model 
  SimulateOpticalModel();
  // Simulate a timing model
  SimulateTimingModel();
  // Merge hits that happened in the same scintillator strip and within the same readout time window
  MergeCoincidentHits();

  EventCounter++;
}

TMS_Event::TMS_Event(TMS_Event &event, int slice) : TMS_Hits(event.GetHits(slice, true)), NonTMS_Hits(event.NonTMS_Hits),
      TMS_TrueParticles(event.TMS_TrueParticles), nTrueForgottenParticles(event.nTrueForgottenParticles),
      TMS_TruePrimaryParticles(event.TMS_TruePrimaryParticles),
      TMS_Tracks(event.TMS_Tracks), Reaction(event.Reaction), Reactions(event.Reactions),
      TrueNeutrino(event.TrueNeutrino), 
      TrueNeutrinoPosition(event.TrueNeutrinoPosition),
      TrueLeptonPosition(event.TrueLeptonPosition), 
      TrueLeptonMomentum(event.TrueLeptonMomentum),  
      TrueVisibleEnergyPerVertex(event.TrueVisibleEnergyPerVertex), 
      TrueVisibleEnergyPerParticle(event.TrueVisibleEnergyPerParticle), 
      ChannelPositions(event.ChannelPositions), 
      DeadChannelTimes(event.DeadChannelTimes), ReadChannelTimes(event.ReadChannelTimes), 
      TimeSliceBounds(event.TimeSliceBounds), info_about_vtx(event.info_about_vtx),
      generator(event.generator) {
  // Create an event from a slice of another event
  SliceNumber = slice;
  SpillNumber = event.SpillNumber;
  
  
  nTrueTrajectories = -999;
  VertexIdOfMostEnergyInEvent = -9991;
  LightWeight = true;
  GetVertexIdOfMostVisibleEnergy();
  
  // Todo, did I copy everything
  // Update event counter if slice != 0, and keep old event number for slice 0.
  if (slice != 0) {
    EventNumber = EventCounter;
    EventCounter++;
  }
  else {
    EventNumber = event.EventNumber;
  }
  
  Reaction = "";
  
  int primary_vertex_id = GetVertexIdOfMostVisibleEnergy();
  if (primary_vertex_id >= 0) {
    SetLeptonInfoUsingVertexID(primary_vertex_id);
    if (Reactions.find(primary_vertex_id) != Reactions.end()) 
      Reaction = Reactions[primary_vertex_id];
    else { Reaction = "NA"; std::cout<<"Warning: couldn't find reaction for primary vertex"<<std::endl; }
  }

  // Update the counts per slice
  ConnectTrueHitWithTrueParticle(true);
}

void TMS_Event::MergeCoincidentHits() {
  std::sort(TMS_Hits.begin(), TMS_Hits.end(), TMS_Hit::SortByZThenT);
  
  const double readout_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();
  
  // Loop over the original hits
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); 
      it != TMS_Hits.end(); it++) {
    if ((*it).GetPedSup()) continue; // Skip hits which are already removed
    // Maybe this hit has already been counted
    double z = (*it).GetZ();
    double y = (*it).GetNotZ();
    //double e = hit.GetE();
    double t = (*it).GetT();

    // Look ahead to find duplicates, but stop when z != z2
    std::vector<std::vector<TMS_Hit>::iterator> duplicates;
    for (std::vector<TMS_Hit>::iterator jt = it+1; jt != TMS_Hits.end() && jt->GetZ() == z; ++jt) {

      TMS_Hit hit2 = *(jt);
      if (hit2.GetPedSup()) continue; // Skip hits which are already removed
      double z2 = hit2.GetZ();
      double y2 = hit2.GetNotZ();
      //double e2 = hit2.GetE();
      double t2 = hit2.GetT();

      // Merge
      if (z == z2 && y == y2 && fabs(t2-t) < readout_time) {
        (*it).MergeWith(hit2);
        // todo, we may want to store an array of true hits. One way would be to move the merging code within the hit class
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
  for (auto& hit : TMS_Hits) {
    if (!hit.GetPedSup()) remaining_hits.push_back(hit);
    else deleted_hits.push_back(hit);
    if (!hit.GetPedSup() && hit.GetE() > 10000)  std::cout << "Warning: Found hit higher than 10 GeV. Seems unlikely. Hit E = " << (hit.GetE() / 1000.0) << " GeV." << std::endl;
  }
  TMS_Hits.clear();
  for (auto& hit : remaining_hits) TMS_Hits.push_back(hit);
  deleted_hits.erase(deleted_hits.begin(), deleted_hits.end());
}

void TMS_Event::SimulateDarkCount() {
  // TODO Add noise hits. They can fake readout
  // One issue is that there's no truth info about the particles to save. 
}

void TMS_Event::SimulateReadoutNoise() {
  // Only want to simulate the little bit of electronic noise from reading out after merging hits
  // Otherwise we're adding together a bunch of random numbers centered around zero, leading to an average of zero
  const double readout_noise = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutNoise();
  if (readout_noise > 0) {
    for (auto& hit : TMS_Hits) {
      double pe = hit.GetPE();
      // Skip 0 pe hits, they'll be removed in ped sup step
      if (pe > 0) {
        double E = hit.GetE();
        // Now take into account electronic readout noise
        std::normal_distribution<double> normal(pe, readout_noise);
        double newpe = normal(generator);
        double newE = E * newpe / pe;
        hit.SetPE(newpe);
        hit.SetE(newE);
      }
    }
  }
}

void TMS_Event::SimulateOpticalModel() {
  // Steps:
  // Loop over hits
  // Convert hit E -> PE
  // Do a poisson throw to get the number of PE
  // Apply some effect of PE capture into the fiber (assumed to be part of E -> PE conversion)
  // Split PE in two randomly for short path vs the long way.

  // TODO add second exponential term using fast decay length
  const double birks_constant = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_BirksConstant(); // mm / MeV
  const bool should_simulate_poisson_throws = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_ShouldSimulatePoisson();
  const bool should_simulate_fiber_lengths = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_ShouldSimulateFiberLengths();
  
  const double wsf_attenuation_length = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFAttenuationLength(); // m
  // In reality, light bounces so there's a length multiplier
  const double wsf_length_multiplier = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFLengthMultiplier();
  const double wsf_decay_constant = 1/wsf_attenuation_length; 
  const double wsf_fiber_reflection_eff = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFEndReflectionEff(); // How much light will reflect at the end
  const double fiber_coupling_eff = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_AdditionalFiberCouplingEff();
  const double optic_fiber_attenuation_length = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_AdditionalFiberAttenuationLength();
  const double optic_fiber_decay_constant = 1/optic_fiber_attenuation_length; 
  const double optic_fiber_length = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_AdditionalFiberLength();
  
  const double readout_coupling_eff = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_ReadoutCouplingEff();
  
  for (auto& hit : TMS_Hits) {
    double pe = hit.GetPE();
    
    // Applies birk's suppression
    double de = hit.GetTrueHit().GetE();
    double dx = hit.GetTrueHit().GetdX();
    double dedx = 0;
    if (dx > 1e-8) dedx = de / dx;
    else dedx = de / 1.0;
    pe *= 1.0 / (1.0 + birks_constant * dedx);
    
    double pe_short = pe;
    double pe_long = 0;
    if (should_simulate_poisson_throws) {
      // Do a poisson throw to get the number of PE
      std::poisson_distribution<int> poisson(pe);
      pe = poisson(generator);
      // Now split the photons into the long and short paths with 50% chance of each
      std::binomial_distribution<int> binomial(pe, 0.5);
      pe_short = binomial(generator);
      pe_long = pe - pe_short;
    }
    
    if (should_simulate_fiber_lengths) {
    
      // Calculate the long and short path lengths
      #ifdef USE_OLD_CODE
      double true_y = hit.GetTrueHit().GetY() / 1000.0; // m
      // In case of orthogonal (X) layers change to GetX()
      if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) true_y = hit.GetTrueHit().GetX() / 1000.0;
      // assuming 0 is center, and assume we're reading out from top, then top would be biased negative and bottom positive, so -true_y.
      // TODO manually found this center. Make function in geom tools that returns values about scint
      // TODO fix math
      double distance_from_middle = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() - true_y;  // -1.54799
      double distance_from_end = distance_from_middle + 2;
      double long_way_distance_from_end = 4 + (4 - distance_from_end);
      #else
      double distance_from_end = hit.GetTrueDistanceFromReadout() * 1e-3; // m
      double long_way_distance_from_end = hit.GetTrueLongDistanceFromReadout() * 1e-3; // m
      #endif
      
      // In reality, light bounces so there's a multiplier
      // TODO it may be more realistic to make this non-linear
      distance_from_end *= wsf_length_multiplier;
      long_way_distance_from_end *= wsf_length_multiplier;
      
      // Now do exponential decay
      pe_short = pe_short * std::exp(-wsf_decay_constant * distance_from_end);
      pe_long = pe_long * std::exp(-wsf_decay_constant * long_way_distance_from_end) * wsf_fiber_reflection_eff;
      
      // Now possibly couple to a regular optical fiber
      if (optic_fiber_length > 0) {
        pe_short = fiber_coupling_eff * pe_short * std::exp(-optic_fiber_decay_constant * optic_fiber_length);
        pe_long = fiber_coupling_eff * pe_long * std::exp(-optic_fiber_decay_constant * optic_fiber_length);
      }
    }
    
    // Now couple between the fibers and the readout
    pe_long = pe_long * readout_coupling_eff;
    pe_short = pe_short * readout_coupling_eff;
    
    // Now save this information
    pe = pe_long + pe_short;
    
    // We want to save info right after fibers but before electronic conversion noise
    // This is particularly useful for timing information which cares about the first photon to be detected
    hit.GetAdjustableTrueHit().SetPEAfterFibers(pe);
    hit.GetAdjustableTrueHit().SetPEAfterFibersLongPath(pe_long);
    hit.GetAdjustableTrueHit().SetPEAfterFibersShortPath(pe_short);

    // Now save the reconstructed information
    hit.SetPE(pe);
    // Need to convert from PE to MeV. Could use 1/LY but have to account for additional effects.
    // so get constant to remove the effect of poisson, birks, fiber length, and readout noise above
    // The largest effect is fiber length
    double calibration_constant = TMS_Manager::GetInstance().Get_RECO_CALIBRATION_EnergyCalibration();
    double reco_e = pe * calibration_constant;  
    hit.SetE(reco_e);  
    hit.SetEVis(reco_e);
  }
}

void TMS_Event::SimulatePedestalSubtraction() {
  // Don't actually remove the hits because they may be relevant for other processes
  // Loop over the hits
  //for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); it++) {
  //   auto hit = (*it);
  for (auto& hit : TMS_Hits) {
     double hit_pe = hit.GetPE();
     if (hit_pe < TMS_Readout_Manager::GetInstance().Get_Sim_Readout_PedestalSubtractionThreshold()) {
       hit.SetPedSup(true);
     }   
  }
}

int TMS_Event::GetUniqIDForDeadtime(const TMS_Hit& hit) const {
  // For a per-channel deadtime, this can return a unique id for a channel
  // But some detectors have deadtime for a whole board. 
  // In that case, this should return a single id for the whole board.
  int id = hit.GetNotZ() + 100000 * hit.GetZ(); // TODO make sure it's unique
  //std::cout<<"x: "<<hit.GetNotZ()<<", z: "<<hit.GetZ()<<", id: "<<id<<std::endl;
  return id;
}

void TMS_Event::SimulateDeadtime() {
  // Simulates readout windows, deadtime and zombie time.
  // |----  readout -----|-------deadtime-------{zombie time}]
  // The actual merging of hits in a readout window happens in MergeCoincidentHits.

  // How long a channel can read out
  double readout_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ReadoutTime();;// TMS_Const::TMS_TimeThreshold; // ns
  // How long a channel is dead before it can read out again
  double deadtime = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_Deadtime();; // ns
  // Imagine a hit before the end of deadtime, but the system is close enough to resetting that you
  // channel starts recording energy. So when readout is ready, it looks like a hit hit at t=readout_ready_time.
  // That's zombie time.
  double zombie_time = TMS_Readout_Manager::GetInstance().Get_Sim_Readout_ZombieTime();;
  
  const bool deadtime_verbose = false;
  
  
  if (deadtime > 0) {
    // Want sorted hits by T so we can find the first hit in a channel to be the start of readout windows and deadtime.
    std::sort(TMS_Hits.begin(), TMS_Hits.end(), TMS_Hit::SortByT);
    
    int n_dead_hits = 0;
    int n_zombie_hits = 0;
    
    // These store the end time for each window by GetUniqIDForDeadtime id.
    // If there's no matching id, then we haven't seen that id yet and this hit can be the start of a readout window.
    std::map<int, double> readout_map;
    std::map<int, double> deadtime_map;
    std::map<int, double> zombie_map;
    std::map<int, bool> has_zombie_map;
    std::map<int, double> x_map;
    std::map<int, double> z_map;
    std::map<int, double> t_map;
    //for (auto& hit : TMS_Hits) {
    for (size_t i = 0; i < TMS_Hits.size(); ++i) {
      auto& hit = TMS_Hits[i];
      double t = hit.GetT();
      const int id = GetUniqIDForDeadtime(hit);
      auto it_read = readout_map.find(id);
      auto it_dead = deadtime_map.find(id);
      auto it_zombie = zombie_map.find(id);
      bool should_reset_times = false;
      if (it_read == readout_map.end()) {
        // This id hasn't been seen before. We can read with no issue
        should_reset_times = true;
        if (deadtime_verbose) std::cout<<"New channel -> Do regular read"<<std::endl;
      }
      else {
        // We have seen this channel. So next we need to see if we're in the readout time or deadtime
        double t_read = it_read->second;
        double t_dead = it_dead->second;
        double t_zombie = it_zombie->second;
        if (deadtime_verbose) std::cout<<"Found channel we found already with Notz: "<<hit.GetNotZ()<<", z: "<<hit.GetZ()<<", id: "<<id<<"\n";
        if (deadtime_verbose) std::cout<<"Compare with previous channel Notz: "<<x_map[id]<<", z: "<<z_map[id]<<", t: "<<t_map[id]<<"\n";
        if (x_map[id] != hit.GetNotZ() || z_map[id] != hit.GetZ()) std::cout<<"\n** Found mismatch in Notz,z **\n"<<std::endl;
        if (deadtime_verbose) std::cout<<"i="<<i<<", t="<<t<<", t_read="<<t_read<<", t_dead="<<t_dead<<", t_zombie="<<t_zombie<<", dt="<<(t-t_read+readout_time)<<", dt_map: "<<(t-t_map[id])<<std::endl;
        if (t < t_read) {
          // We can do a regular read
          if (deadtime_verbose) std::cout<<"t < t_read -> Do regular read"<<std::endl;
        } 
        else if (zombie_time > 0 && t < t_zombie) {
          // Zombie time sets the time to the start of the next readout which is the end of deadtime
          hit.SetT(t_dead);
          n_zombie_hits += 1;
          has_zombie_map[id] = true;
          if (deadtime_verbose) std::cout<<"t < t_zombie -> Treat as zombie"<<std::endl;
        }
        else if (t < t_dead) {
          // Suppress channel since it's dead
          hit.SetPedSup(true);
          n_dead_hits += 1;
          if (deadtime_verbose) std::cout<<"t < t_dead -> Kill channel"<<std::endl;
        }
        else {
          // This channel is past the deadtime so reset our windows
          if (deadtime_verbose) std::cout<<"t >= t_dead -> Regular read and reset windows"<<std::endl;
          should_reset_times = true;
        }
      }
      
      // Calculate the times that this id can be read or is dead or is zombie
      if (should_reset_times) {
        // If there's a zombie time, that hit should be the start of the next readout
        // And its time will be the end of the current deadtime
        // But we don't need to do this check if this hit is outside the deadtime window of the next hit
        // So calculate that window first
        // But also this hit now needs to be checked against the zombie times's windows, so redo this hit
        double deadtime_window_starting_from_end_of_deadtime = it_dead->second + deadtime;
        if (deadtime_map.find(id) != deadtime_map.end() && has_zombie_map[id] == true && 
          t < deadtime_window_starting_from_end_of_deadtime) { 
          t = deadtime_map[id];
          has_zombie_map[id] = false;
          // Need to redo this hit to check that it isn't in the deadtime of the zombie hit
          i--;
        }
        // Since this is the first hit for an id, the end of the readout window is t + readout_time
        double t_read = t + readout_time;
        readout_map[id] = t_read;
        // Deadtime starts at the end of the readout period
        double t_dead = t_read + deadtime;
        deadtime_map[id] = t_dead;
        // Zombie time is the time before the end of deadtime
        // So 20ns of zombie time in a 500ns deadtime would start at 480ns.
        double t_zombie = t_dead - zombie_time;
        zombie_map[id] = t_zombie;
        
        x_map[id] = hit.GetNotZ();
        z_map[id] = hit.GetZ();
        t_map[id] = t;
        
        auto position = std::make_pair(hit.GetNotZ(), hit.GetZ());
        auto deadtime_range = std::make_pair(t_read, t_dead);
        auto readout_range = std::make_pair(t, t_read);
        ChannelPositions.push_back(position);
        DeadChannelTimes.push_back(deadtime_range);
        ReadChannelTimes.push_back(readout_range);
        
        #ifdef RECORD_HIT_DEADTIME
        // In theory merging hits should capture correct deadtimes based on the first deadtime recorded here
        // But if doing things by sipm or some larger group, then we're not recording all the info
        hit.SetDeadtimeStart(t_read);
        hit.SetDeadtimeStop(t_dead);
        #endif
      } 
      
      // TODO remove
      //if (i > 100) exit(0);
      
    }
    double n_dead_hits_as_percent = 100.0 * n_dead_hits / (float)TMS_Hits.size();
    std::cout<<"N dead hits: "<<n_dead_hits<<" out of "<<TMS_Hits.size()<<" hits. That's "<<n_dead_hits_as_percent<<"%"<<std::endl;
    if (zombie_time > 0) std::cout<<"N zombie hits: "<<n_zombie_hits<<std::endl;
  }
}

void TMS_Event::SimulateTimingModel() {
  // List of timing effects to simulate:
  // Random electronic timing noise
  // Deadtime
  // Time skew from first PE to hit sensor
  // Optical fiber length delays (corrected to strip center)
  // Timing effects from random noise, cross talk, afterpulsing
  // TODO check constants or put in config  
  std::normal_distribution<double> noise_distribution(0.0, 1); // Mean of 0.0 and standard deviation of 1ns
  std::uniform_int_distribution<int> coin_flip(0, 1); // 0 or 1 depending on if you went long or short
  double scintillator_decay_time = 3.0; // ns
  double wsf_decay_time = 20.0; // ns
  std::exponential_distribution<double> exp_scint(1 / scintillator_decay_time); // Decay time = 3ns for scintillator
  std::exponential_distribution<double> exp_wsf(1 / wsf_decay_time); // 20ns for wavelength shifting fiber
  const double SPEED_OF_LIGHT =  0.2998; // m/ns
  const double FIBER_N = 1.5; // 
  const double SPEED_OF_LIGHT_IN_FIBER = SPEED_OF_LIGHT / FIBER_N;
  
  const double wsf_length_multiplier = TMS_Readout_Manager::GetInstance().Get_Sim_Optical_WSFLengthMultiplier();
  
  //double avg = 0;
  //double maxy = -1e9;
  //double miny = 1e9;
  //int n = 0;
  for (auto& hit : TMS_Hits) {
    double t = 0;
    // Random electronic timing noise (~1ns or less)
    t += noise_distribution(generator);
    // Optical fiber length delay (corrected to strip center) 
    // (up to 13.4ns assuming 4m from edge, but correlated with y position. If delta y = 1m spread, than relative error is only 3.3ns)
    #ifdef USE_OLD_CODE
    double true_y = hit.GetTrueHit().GetY() / 1000.0; // m
    // Making sure this gets changed for orthogonal (X) layers
    if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) true_y = hit.GetTrueHit().GetX() / 1000.0;
    //miny = std::min(miny, true_y);
    //maxy = std::max(maxy, true_y);
    // assuming 0 is center, and assume we're reading out from top, then top would be biased negative and bottom positive, so -true_y.
    // TODO manually found this center. Want a better way in case things change
    double distance_from_middle = TMS_Manager::GetInstance().Get_Geometry_YMIDDLE() - true_y;  //-1.54799 
    double long_way_distance = distance_from_middle + 8;
    #else
    double distance_from_middle = hit.GetTrueDistanceFromMiddle() * 1e-3; // m
    double long_way_distance = hit.GetTrueLongDistanceFromMiddle() * 1e-3; // m
    #endif
    
    // In reality, light bounces so there's a multiplier to the distance
    // todo, it may be more realistic to make this non-linear
    distance_from_middle *= wsf_length_multiplier;
    long_way_distance *= wsf_length_multiplier;
    
    // Find the time correction
    double time_correction = distance_from_middle / SPEED_OF_LIGHT_IN_FIBER;
    // This is the time correction if you go the long way instead
    double time_correction_long_way = long_way_distance / SPEED_OF_LIGHT_IN_FIBER;
    
    // Time slew (up to 30ns for 1pe hits, 9ns for 5pe, ~2ns 22pe. Typically 22pe mips assuming 45 pe mips with half going the long way)
    double pe_short_path = hit.GetTrueHit().GetPEAfterFibersShortPath();
    double pe_long_path = hit.GetTrueHit().GetPEAfterFibersLongPath();
    double minimum_time_offset = 1e100;
    
    #define USE_GAMMA_DISTRIBUTION
    #ifdef USE_GAMMA_DISTRIBUTION
    // TODO test this version with gamma
    
    // Gamma distribution does throws without having to do the throws
    // So it answer the question, what's the lowest of a random exponential assuming N throws.
    std::gamma_distribution<double> gamma_scint_short_path(pe_short_path, 1.0 / scintillator_decay_time);
    std::gamma_distribution<double> gamma_wsf_short_path(pe_short_path, 1.0 / wsf_decay_time);
    double minimum_time_gamma_scint_short_path = gamma_scint_short_path(generator);
    double minimum_time_gamma_wsf_short_path = gamma_wsf_short_path(generator);
    double minimum_time_offset_short_path = minimum_time_gamma_scint_short_path + minimum_time_gamma_wsf_short_path + time_correction;
    minimum_time_offset = std::min(minimum_time_offset_short_path, minimum_time_offset);
    
    std::gamma_distribution<double> gamma_scint_long_path(pe_long_path, 1.0 / scintillator_decay_time);
    std::gamma_distribution<double> gamma_wsf_long_path(pe_long_path, 1.0 / wsf_decay_time);
    double minimum_time_gamma_scint_long_path = gamma_scint_long_path(generator);
    double minimum_time_gamma_wsf_long_path = gamma_wsf_long_path(generator);
    double minimum_time_offset_long_path = minimum_time_gamma_scint_long_path + minimum_time_gamma_wsf_long_path + time_correction_long_way;
    minimum_time_offset = std::min(minimum_time_offset_long_path, minimum_time_offset);
    
    #elif
    // We don't have to do 1000s of throws. The time will be very close to zero.
    // Assuming 1k PE, the mean time is ~0.02ns vs ~0.06ns for 300 PE.
    const double MAX_PE_THROWS = 300;
    if (pe_short_path > MAX_PE_THROWS) {
      pe_short_path = MAX_PE_THROWS;
    }
    while (pe_short_path > 0) {
      double time_offset = time_correction;
      time_offset += exp_scint(generator);
      time_offset += exp_wsf(generator);
      minimum_time_offset = std::min(time_offset, minimum_time_offset);
      pe_short_path -= 1;
    }
    if (pe_long_path > MAX_PE_THROWS) {
      pe_long_path = MAX_PE_THROWS;
    }
    while (pe_long_path > 0) {
      double time_offset = time_correction_long_way;
      time_offset += exp_scint(generator);
      time_offset += exp_wsf(generator);
      minimum_time_offset = std::min(time_offset, minimum_time_offset);
      pe_long_path -= 1;
    }
    #endif
    t += minimum_time_offset;
    double hit_time = hit.GetT();
    //std::cout<<"dt: "<<t<<", hit t: "<<hit_time<<", reco t: "<<hit_time + t<<", min t offset: "<<minimum_time_offset<<", t corr: "<<time_correction<<", dist from middle: "<<distance_from_middle<<", long way t corr: "<<time_correction_long_way<<", long way dist: "<<long_way_distance<<", hit pe: "<<hit.GetPE()<<std::endl;
    //std::cout<<"Hit time: "<<hit_time<<std::endl;
    //std::cout<<"Adjusted hit time: "<<hit_time + t<<std::endl;
    // Finally set the time
    hit.SetT(hit_time + t);
  }
  //avg /= n;
  //std::cout<<"Avg middle: "<<avg<<std::endl;
  /*std::cout<<"Max y: "<<maxy<<std::endl;
  std::cout<<"Min y: "<<miny<<std::endl;
  std::cout<<"Center y: "<<0.5*(maxy - miny)<<std::endl;*/
}

void TMS_Event::ApplyReconstructionEffects() {
  // First apply energy and timing models. Then merge hits. Then do a pedestal subtraction.
  // Simulate an optical model 
  //SimulateOpticalModel();
  // Noise hits can simulate hits
  SimulateDarkCount();
  // Simulate a timing model
  //SimulateTimingModel();
  // Simulate deadtime if needed
  SimulateDeadtime();
  // Merge hits that happened in the same scintillator strip and within the same readout time window
  MergeCoincidentHits();
  // After merging hits, we have a single readout. This readout will have some electronic noise
  SimulateReadoutNoise();
  // Simulate pedestal subtraction where any hit under Get_Sim_Readout_PedestalSubtractionThreshold is removed
  SimulatePedestalSubtraction();
}

const std::vector<TMS_Hit> TMS_Event::GetHits(int slice, bool include_ped_sup) {
  std::vector<TMS_Hit> out;
  for (auto& hit : TMS_Hits) {
    if (!hit.GetPedSup() || include_ped_sup) {
      int slice_number = hit.GetSlice();
      //std::cout<<"Slice number for hit: "<<slice_number<<std::endl;
      if (slice_number == slice || slice < 0) out.push_back(hit);
    }
  }
  return out;
}

bool TMS_Event::IsInTimeSlice(double time) const {
  int current_time_slice = GetSliceNumber();
  bool out;
  if (current_time_slice == 0) {
    // Special case: Make sure t isn't part of any other time slice
    out = true;
    for (const auto& bounds : TimeSliceBounds) {
      double start = bounds.first;
      double end = bounds.second;
      // If t is within any bound, then it's not part of slice zero so it's not in slice 0
      if (start <= time && time <= end) { out = false; break; }
    }
  }
  else {
    // Check if t is within time slice bounds
    if (current_time_slice < 0 || current_time_slice > (int) TimeSliceBounds.size()) {
      std::cout<<"Fatal: IsInTimeSlice got slice number outside time slice bounds. Got: "<<current_time_slice;
      std::cout<<", TimeSliceBounds.size(): "<<TimeSliceBounds.size()<<std::endl;
      throw std::runtime_error("Fatal: IsInTimeSlice got slice number outside time slice bounds");
    }
    double start = TimeSliceBounds[current_time_slice].first;
    double end = TimeSliceBounds[current_time_slice].second;
    if (start <= time && time <= end) out = true;
    else out = false;
  }
  return out;
}

std::pair<double, double> TMS_Event::GetTimeSliceBounds(int slice) {
  if (slice == -1) {
    if (GetSliceNumber() == -1) {
      std::cout<<"Warning: Found circular logic in GetTimeSliceBounds. Returning default bounds"<<std::endl;
      return std::make_pair(0.0, 10000.0);
    }
    return GetTimeSliceBounds(GetSliceNumber());
  }
  else {
    if (slice < 0 || slice >= (int) TimeSliceBounds.size()) {
      std::cout<<"Fatal: GetTimeSliceBounds error: slice: "<<slice<<" outside TimeSliceBounds.size(): "<<TimeSliceBounds.size()<<std::endl;
      throw std::runtime_error("GetTimeSliceBounds error: slice outside TimeSliceBounds range");
    }
    return TimeSliceBounds[slice];
  }
}

// Add a separate event to this event
// Handy for making hacked overlays
void TMS_Event::AddEvent(TMS_Event &Other_Event) {
  //std::cout << "Adding event " << Other_Event.GetEventNumber() << " to event " << GetEventNumber() << std::endl;

  // Get the other hits
  std::vector<TMS_Hit> other_hits = Other_Event.GetHits(-1, true);

  // And use them to expand on the original hits in the event
  for (auto &hit: other_hits) {
    TMS_Hits.emplace_back(std::move(hit));
  }
  
  // Do the same for non-tms hits
  for (auto &hit: Other_Event.NonTMS_Hits) {
    NonTMS_Hits.emplace_back(std::move(hit));
  }

  // Do the same for the true particles
  std::vector<TMS_TrueParticle> other_truepart = Other_Event.GetTrueParticles();
  for (auto &part: other_truepart) {
    TMS_TrueParticles.emplace_back(std::move(part));
  }
  // And true primary particles
  std::vector<TMS_TrueParticle> other_trueprimpart = Other_Event.TMS_TruePrimaryParticles;
  for (auto &part: other_trueprimpart) {
    TMS_TruePrimaryParticles.emplace_back(std::move(part));
  }
  
  // Merge these lists
  TrueVisibleEnergyPerVertex.merge(Other_Event.TrueVisibleEnergyPerVertex);
  TrueVisibleEnergyPerParticle.merge(Other_Event.TrueVisibleEnergyPerParticle);
  Reactions.merge(Other_Event.Reactions);
  info_about_vtx.merge(Other_Event.info_about_vtx);
  // Reset this to recalculate on next call
  VertexIdOfMostEnergyInEvent = -9999;

  nVertices += Other_Event.nVertices;
}

void TMS_Event::OverlayEvents(std::vector<TMS_Event>& overlay_events) {
  for (auto &event : overlay_events) AddEvent(event);
}

void TMS_Event::FinalizeEvent() {
  // Apply the det sim now, after overlaying events
  // The timing and optical model were moved to the initial event creation
  ApplyReconstructionEffects();
  // Connect true vis E and true n hits with true particles
  ConnectTrueHitWithTrueParticle(false);
}

// For now just fill the true neutrino
// But shows how you can easily make a vector of rootracker particles for the TMS_Event to carry along
void TMS_Event::FillTruthFromGRooTracker(int pdg[__EDEP_SIM_MAX_PART__], double p4[__EDEP_SIM_MAX_PART__][4], 
  double vtx[__EDEP_SIM_MAX_PART__][4]) {
  // Momenta/Energy are in GeV
  TrueNeutrino.first.SetX(p4[0][0]);
  TrueNeutrino.first.SetY(p4[0][1]);
  TrueNeutrino.first.SetZ(p4[0][2]);
  TrueNeutrino.first.SetT(p4[0][3]);
  TrueNeutrino.second = pdg[0];
  // Positions are in m
  TrueNeutrinoPosition.SetX(vtx[0][0]);
  TrueNeutrinoPosition.SetY(vtx[0][1]);
  TrueNeutrinoPosition.SetZ(vtx[0][2]);
  TrueNeutrinoPosition.SetT(vtx[0][3]);
  
  if (info_about_vtx.size() == 1) {
    auto it = info_about_vtx.begin();
    auto key = (*it).first;
    auto second = (*it).second;
    // Calculate the distance squared to make sure they're about the same vertex position
    double mm = 1000.0; // convert from m
    double dist2 = 0;
    dist2 += pow(vtx[0][0]*mm - second.vtx.X(), 2);
    dist2 += pow(vtx[0][1]*mm - second.vtx.Y(), 2);
    dist2 += pow(vtx[0][2]*mm - second.vtx.Z(), 2);
    double eps = 0.1; // should be about the same
    if (dist2 < eps) {
      info_about_vtx[key].pdg = pdg[0];
      double MeV = 1000.0; // Convert from GeV
      info_about_vtx[key].p4 = TLorentzVector(p4[0][0] * MeV, p4[0][1] * MeV, p4[0][2] * MeV, p4[0][3] * MeV);
    }
    else { 
      std::cout<<"Found mismatch between groo vtx and info_about_vtx. Please investigate. Dist2="<<dist2<<std::endl;
      std::cout<<vtx[0][0]<<", "<<vtx[0][1]<<", "<<vtx[0][2]<<std::endl;
      std::cout<<second.vtx.X()<<", "<<second.vtx.Y()<<", "<<second.vtx.Z()<<std::endl;
    }
  }
}

void TMS_Event::FillTrueLeptonInfo(int pdg, TLorentzVector position, TLorentzVector momentum, int vertexid) {
  TrueLeptonPDG = pdg;
  TrueLeptonPosition = position;
  TrueLeptonMomentum = momentum;
  TrueLeptonVertexID = vertexid;
}

int TMS_Event::GetVertexIdOfMostVisibleEnergy() {
  // Return early if we've already calculated it
  if (VertexIdOfMostEnergyInEvent >= 0) return VertexIdOfMostEnergyInEvent;

  // Reset the map
  TrueVisibleEnergyPerVertex.clear();
  int min_vertex_id_seen = 1e9;
  int max_vertex_id_seen = -1;
  // First find how much energy is in each variable
  for (auto& hit : TMS_Hits) {
    for (size_t i = 0; i < hit.GetTrueHit().GetNTrueParticles(); i++) {
      int vertex_id = hit.GetTrueHit().GetVertexIds(i);
      // todo, true or reco energy?
      double energy = hit.GetTrueHit().GetEnergyShare(i);
      //std::cout<<"vertex_id = "<<vertex_id<<", energy = "<<energy<<std::endl;
      TrueVisibleEnergyPerVertex[vertex_id] += energy;
      if (vertex_id < min_vertex_id_seen) min_vertex_id_seen = vertex_id;
      if (vertex_id > max_vertex_id_seen) max_vertex_id_seen = vertex_id;
    }
  }
  
  // Now find the most energetic vertex
  double max = -1e9;
  int max_vertex_id = -9992;
  double total_energy = 0;
  for (auto it : TrueVisibleEnergyPerVertex) {
    double vertex = it.first;
    double energy = it.second;
    if (energy >= max) { max = energy; max_vertex_id = vertex; }
    total_energy += energy;
  }
  VertexIdOfMostEnergyInEvent = max_vertex_id;
  VisibleEnergyFromVertexInSlice = max;
  VisibleEnergyFromOtherVerticesInSlice = total_energy - max;
  
  if (TrueVisibleEnergyPerVertex.find(VertexIdOfMostEnergyInEvent) != TrueVisibleEnergyPerVertex.end())
    TotalVisibleEnergyFromVertex = TrueVisibleEnergyPerVertex[VertexIdOfMostEnergyInEvent];
  else {
    if (TrueVisibleEnergyPerVertex.size() > 0) {
      std::cout<<"Warning in GetVertexIdOfMostVisibleEnergy: TrueVisibleEnergyPerVertex.Find(VertexIdOfMostEnergyInEvent) == TrueVisibleEnergyPerVertex.end()"<<std::endl;
    }
  }
  
  return VertexIdOfMostEnergyInEvent;
}

std::pair<double, double> TMS_Event::GetEventTimeRange() {
  double min_time = 1e9;
  double max_time = -1e9;
  for (auto& hit : TMS_Hits) {
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
    for (auto pnt = pos.begin(); (pnt+1) != pos.end(); ++pnt, ++num) {
      auto nextpnt = *(pnt+1);
      TVector3 point1((*pnt).X(), (*pnt).Y(), (*pnt).Z());  //-200
      TVector3 point2(nextpnt.X(), nextpnt.Y(), nextpnt.Z()); //-200
      if (TMS_Geom::GetInstance().IsInsideTMS(point1) && TMS_Geom::GetInstance().IsInsideTMS(point2)) {
        if ((point2-point1).Mag() > 100) {
          continue;
        }
        double tracklength = TMS_Geom::GetInstance().GetTrackLength(point1, point2);
        total += tracklength;
      }
    }
  }
  return total;
}

int TMS_Event::GetTrueParticleIndex(int vertexid, int trackid) {
  int out = -1;
  // Gracefully deal with trackid = -999
  if (vertexid >= 0 && trackid >= 0) {
    for (size_t i = 0; i < TMS_TrueParticles.size(); i++) {
      auto& tp = TMS_TrueParticles.at(i);
      if (tp.GetVertexID() == vertexid && tp.GetTrackId() == trackid) { out = i; break; }
    }
  }
  else {
    std::cout<<"GetTrueParticleIndex: Case of vertex < 0 or trackid < 0. Vertex id: "<<vertexid<<", track id: "<<trackid<<", n TMS_TrueParticles: "<<TMS_TrueParticles.size()<<std::endl;
  }
  if (out < 0) std::cout<<"GetTrueParticleIndex: Case where out < 0. Vertex id: "<<vertexid<<", track id: "<<trackid<<", n TMS_TrueParticles: "<<TMS_TrueParticles.size()<<std::endl;
  return out;
}

int TMS_Event::GetPrimaryLeptonOfVertexID(int vertexid) {
  int lepton_index = -999;
  int current_index = 0;
  for (auto& particle : TMS_TruePrimaryParticles) {
    if (particle.GetVertexID() == vertexid) {
      int pdg = std::abs(particle.GetPDG());
      if (pdg >= 11 && pdg <= 16) {
        lepton_index = current_index;
        break;
      }
    }
    current_index += 1;
  }
  return lepton_index;
}

void TMS_Event::SetLeptonInfoUsingVertexID(int vertexid) {

  // And now fill lepton info
  auto particle_index = GetPrimaryLeptonOfVertexID(vertexid);
  if (particle_index >= 0) {
    auto particle = TMS_TruePrimaryParticles[particle_index];
    int lepton_pdg = particle.GetPDG();
    auto lepton_position = particle.GetBirthPosition();
    auto lepton_momentum = particle.GetBirthMomentumAsLorentz();
    FillTrueLeptonInfo(lepton_pdg, lepton_position, lepton_momentum, vertexid);
  }
  else {
    std::cout<<"Warning in SetLeptonInfoUsingVertexID: GetPrimaryLeptonOfVertexID didn't"
               "return a valid particle index for vertex id "<<vertexid<<std::endl;
    FillTrueLeptonInfo(-9999999, TLorentzVector(-9999999, -999999, -999999, -999999), 
      TLorentzVector(-9999999, -999999, -999999, -999999), vertexid);
  }
}

double TMS_Event::CalculateEnergyInLArOuterShell(double thickness, int vertexid) {
  double out = 0;
  // Lar doesn't have good timing info, so we want all non tms hits, not just in this slice
  for (const auto& hit : NonTMS_Hits) {
    if (vertexid < 0 || hit.GetVertexIds(0) == vertexid) {
      TVector3 position(hit.GetX(), hit.GetY(), hit.GetZ());
      if (TMS_Geom::GetInstance().IsInsideLAr(position) && !TMS_Geom::GetInstance().IsInsideLAr(position, thickness)) {
        out += hit.GetHadronicEnergy();
      }
    }
  }
  return out;
}

double TMS_Event::CalculateEnergyInLAr(int vertexid) {
  double out = 0;
  for (const auto& hit : NonTMS_Hits) {
    if (hit.GetVertexIds(0) < 0) std::cout<<"Warning: found true hit with < 0 VertexId"<<std::endl;
    if (vertexid < 0 || hit.GetVertexIds(0) == vertexid) { 
      TVector3 position(hit.GetX(), hit.GetY(), hit.GetZ());
      if (TMS_Geom::GetInstance().IsInsideLAr(position))
        out += hit.GetHadronicEnergy();
    }
  }
  return out;
}


double TMS_Event::CalculateTotalNonTMSEnergy(int vertexid) {
  double out = 0;
  for (const auto& hit : NonTMS_Hits) {
    if (hit.GetVertexIds(0) < 0) std::cout<<"Warning: found true hit with < 0 VertexId"<<std::endl;
    if (vertexid < 0 || hit.GetVertexIds(0) == vertexid) out += hit.GetHadronicEnergy();
  }
  return out;
}

void TMS_Event::ConnectTrueHitWithTrueParticle(bool slice) {
  // Now count the number of true hits per particle
  std::map<int, int> NHitsPerParticle;
  std::map<int, double> EnergyPerParticle;
  for (auto& hit : TMS_Hits) {
    // Only count hits that are not ped subtracted
    if (!hit.GetPedSup()) {
      auto true_hit = hit.GetTrueHit();
      // Only add 1 hit for each key once, so track if we saw a key already
      std::map<int, int> key_seen;
      for (size_t i = 0; i < true_hit.GetNTrueParticles(); i++) {
        int key = true_hit.GetVertexIds(i) * 100000 + true_hit.GetPrimaryIds(i);
        if (key_seen.find(key) == key_seen.end()) { 
          NHitsPerParticle[key] += 1;
          key_seen[key] = 1;
        }
        EnergyPerParticle[key] += true_hit.GetEnergyShare(i);
      }
    }
  }
  for (size_t i = 0; i < TMS_TrueParticles.size(); i++) {
    int count = 0;
    double energy = 0;
    // If it's not in the map, don't create it
    int key = TMS_TrueParticles[i].GetVertexID() * 100000 + TMS_TrueParticles[i].GetTrackId();
    if (NHitsPerParticle.find(key) != NHitsPerParticle.end()) {
      count = NHitsPerParticle[key];
      energy = EnergyPerParticle[key];
    }
    TMS_TrueParticles[i].SetNTrueHits(count, slice);
    TMS_TrueParticles[i].SetTrueVisibleEnergy(energy, slice);
  }
}


void TMS_Event::SaveKeyVertexInfo(const TMS_TrueHit& hit) {
  for (size_t i = 0; i < hit.GetNTrueParticles(); i++) {
    int vertex_id = hit.GetVertexIds(i);
    if (info_about_vtx.find(vertex_id) != info_about_vtx.end()) {
      info_about_vtx[vertex_id].AddEnergyFromHit(hit, i);
    }
    else std::cout<<"This should not happen but I didn't find a vertex for vertex id "<<vertex_id<<std::endl;
  }
}

Vtx_Info* TMS_Event::GetVertexInfo(int vertex_id) {
  Vtx_Info* out = NULL;
  if (info_about_vtx.find(vertex_id) != info_about_vtx.end())
    out = &info_about_vtx.at(vertex_id);
  return out;
}


void Vtx_Info::AddEnergyFromHit(const TMS_TrueHit& hit, int index) {
  double hadronic_energy = hit.GetHadronicEnergy() * hit.GetEnergySharePortion(index);
  double energy = hit.GetE() * hit.GetEnergySharePortion(index);
  TVector3 position(hit.GetX(), hit.GetY(), hit.GetZ());

  // Total
  hadronic_energy_total += hadronic_energy;
  true_visible_energy_total += energy;

  // Lar-specific
  if (TMS_Geom::GetInstance().IsInsideLAr(position)) {
    hadronic_energy_lar += hadronic_energy;
    true_visible_energy_lar += energy;
  }
  // Lar outer-shell for the hadron containment cut
  if (TMS_Geom::GetInstance().IsInsideLArShell(position) && 0 < hadronic_energy) {
    hadronic_energy_lar_shell += hadronic_energy;
    UpdateShellEnergyCut();
  }

  // TMS-specific
  if (TMS_Geom::GetInstance().IsInsideLAr(position)) {
    hadronic_energy_tms += hadronic_energy;
    true_visible_energy_tms += energy;
  }
}






