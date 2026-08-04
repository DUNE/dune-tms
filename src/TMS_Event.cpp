#include "TMS_Event.h"
#include "TMS_Readout_Manager.h"
#include "TMS_VertexId.h"
#include "TMS_DetectorSimulation.h"
#include "TMS_SignalProcessing.h"
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
  VertexGlobalIdOfMostEnergyInEvent = -999;
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
  if (event.Primaries.size() > 1) {
    std::cout<<"Fatal: TMS_Event found "<<event.Primaries.size()
             <<" primary vertices in one TG4Event. This path used to fall back to"
             <<" PrimaryVertex::GetInteractionNumber(), but that is not guaranteed to match"
             <<" the edep-sim EventId convention used by CAFMaker."<<std::endl;
    throw std::runtime_error("Fatal: multiple primary vertices need explicit edep-sim vertex IDs");
  }
  // Loop over the primary vertices
  for (TG4PrimaryVertexContainer::iterator it = event.Primaries.begin(); it != event.Primaries.end(); ++it) {

    TG4PrimaryVertex vtx = *it;
    Reaction = (*it).GetReaction();
    
    // Interaction number is off-by-one in recent microprod files, so set it manually
    // See https://github.com/DUNE/2x2_sim/issues/61
    vtx.InteractionNumber = current_vertexid;
    if (current_vertexid < 0) {
      std::cout<<"Fatal: Got a current_vertexid < 0 in TMS_Event: "<<current_vertexid<<std::endl;
      throw std::runtime_error("Fatal: Get a vertex id < 0");
    }
      
    Vtx_Info vtx_info;
    vtx_info.reaction = Reaction;
    vtx_info.run_id = event.RunId;
    vtx_info.vtx_id = current_vertexid;
    // Had issues with lorentz vectors before so best make a copy
    vtx_info.SetVtx(TLorentzVector(vtx.GetPosition().X(), vtx.GetPosition().Y(), vtx.GetPosition().Z(), vtx.GetPosition().T()));
    
    const long long global_vertex_id = TMS_MakeGlobalVertexID(vtx_info.run_id, vtx_info.vtx_id);
    Reactions[global_vertex_id] = Reaction;
    auto inserted = info_about_vtx.emplace(global_vertex_id, vtx_info);
    if (!inserted.second) {
      std::cout<<"Fatal: Duplicate global vertex id "<<global_vertex_id
               <<" for run "<<vtx_info.run_id<<" vertex "<<vtx_info.vtx_id<<std::endl;
      throw std::runtime_error("Fatal: duplicate global vertex id");
    }

    if (FillEvent) {
      // Primary particles in edep-sim are before any particle propagation happens
      // i.e. it's the particles out of the neutrino generation; don't save them
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;

      // Loop over the particles in the vertex and save them
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
        TMS_TrueParticle truepart = TMS_TrueParticle(particle, vtx, event.RunId);
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
            TMS_Manager &manager = TMS_Manager::GetInstance();
            // Check the TMS volume first
            if (VolumeName.find(manager.Get_GEOMETRY_VOLUME_TMSVolume()) == std::string::npos &&
                VolumeName.find(manager.Get_GEOMETRY_VOLUME_ModuleLayer()) == std::string::npos &&
                VolumeName.find(manager.Get_GEOMETRY_VOLUME_TMSEDepSimVolume()) == std::string::npos) continue;

            // check the LAr volume
            if (TMSLArOnly) {
              if (VolumeName.find(manager.Get_GEOMETRY_VOLUME_LArActive()) == std::string::npos) continue;
            }
          }

          // If firsttime is true and the above passes, this is the first time the particle enters any volume of interest, so create it
          if (firsttime) {
            // Can't set start momentum and position whe looping over the trajectory points, do this later
            //TMS_TrueParticle part(ParentId, TrackId, PDGcode, Momentum, Position);
            TMS_TrueParticle part(ParentId, TrackId, PDGcode, current_vertexid, event.RunId);
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
  std::map<int, long long> mapping_track_to_vertex_global_id;
  const int vertex_index = event.EventId;
  const long long vertex_global_index = TMS_MakeGlobalVertexID(RunNumber, vertex_index);
  for (auto vertex : event.Primaries) {
    for (auto particle : vertex.Particles) {
      int track_id = particle.GetTrackId();
      mapping_track_to_vertex_global_id[track_id] = vertex_global_index;
    }
    for (auto traj : event.Trajectories) {
      int track_id = traj.GetTrackId();
      mapping_track_to_vertex_global_id[track_id] = vertex_global_index;
    }
  }

  std::map<std::pair<long long, int>, TMS_TrueParticle*> mapping_track_to_true_particle;
  for (auto& tp : TMS_TrueParticles) {
    auto key = std::make_pair(TMS_MakeGlobalVertexID(tp.GetRunID(), tp.GetVertexID()), tp.GetTrackId());
    mapping_track_to_true_particle[key] = &tp;
  }
  std::map<std::tuple<int, int, int, long long>, size_t> map_pos_nontms_hits;

  // Loop over each hit
  for (TG4HitSegmentDetectors::iterator jt = event.SegmentDetectors.begin(); jt != event.SegmentDetectors.end(); ++jt) {
    // Only look at TMS hits
    std::string DetString = (*jt).first;

    // Skip hits outside of the TMS if running lightweight
    if (TMSOnly && DetString != TMS_Manager::GetInstance().Get_GEOMETRY_VOLUME_TMSEDepSimVolume()) continue;

    TG4HitSegmentContainer tms_hits = (*jt).second;
    for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
      TG4HitSegment edep_hit = *kt;
      int track_id = edep_hit.GetPrimaryId();
      long long vertex_global_id = -999;
      auto value = mapping_track_to_vertex_global_id.find(track_id);
      if (value == mapping_track_to_vertex_global_id.end()) {
        std::cout<<"WARNING: Didn't find track id in mapping_track_to_vertex_global_id! track_id = "<<track_id<<", mapping_track_to_vertex_global_id.size() = "<<mapping_track_to_vertex_global_id.size()<<", this shouldn't happen anymore\n\n\n"<<std::endl;
      }
      else vertex_global_id = value->second;
      TMS_Hit hit = TMS_Hit(edep_hit);
      hit.SetHitId(NextHitId());
      int barnum = hit.GetBarNumber();
      // Only add if within the TMS
      // Can't use x,y or z because geometry might change. But we know things aren't set if there's no bar number
      if (barnum >= 0) {
        // Truth is constructed separately from the reco-level TMS_Hit above (Phase III --
        // TMS_TrueHit is no longer embedded in TMS_Hit) and stored in the event-level side
        // table, keyed by this hit's HitId.
        TMS_TrueHit t(edep_hit, vertex_global_id);
        for (size_t i = 0; i < t.GetNTrueParticles(); i++) {
          auto key = std::make_pair(t.GetVertexGlobalIds(i), t.GetPrimaryIds(i));
          if (mapping_track_to_true_particle.find(key) != mapping_track_to_true_particle.end()) {
            // Now set info
            auto tp = mapping_track_to_true_particle[key];
            if (tp->IsLeptonic()) t.SetEnergyLeptonic(i);
          }
        }
        SaveKeyVertexInfo(t);
        SetTrueHit(hit.GetHitId(), t);
        TMS_Hits.push_back(std::move(hit));

        // todo, maybe skip for michel electrons or late neutrons
        for (size_t i = 0; i < t.GetNTrueParticles(); i++) {
          TrueVisibleEnergyPerVertex[t.GetVertexGlobalIds(i)] += t.GetEnergyShare(i);
          TrueVisibleEnergyPerParticle[std::make_pair(t.GetVertexGlobalIds(i), t.GetPrimaryIds(i))] += t.GetEnergyShare(i);
        }
      }
      else if (DetString.find(TMS_Manager::GetInstance().Get_GEOMETRY_VOLUME_LArActive()) != std::string::npos) {
        // Only care about LAr active volume
        // We only need it for truth info so just save truth info
        TMS_TrueHit t(edep_hit, vertex_global_id);
        for (size_t i = 0; i < t.GetNTrueParticles(); i++) {
          auto key = std::make_pair(t.GetVertexGlobalIds(i), t.GetPrimaryIds(i));
          auto itp = mapping_track_to_true_particle.find(key);
          if (itp != mapping_track_to_true_particle.end()) {
            // Now set info
            auto tp = itp->second;
            if (tp->IsLeptonic()) t.SetEnergyLeptonic(i);
          }
        }
        
        SaveKeyVertexInfo(t);
        
        double divide = 10.0;
        auto poskey = std::tuple((int) (t.GetX() / divide), (int) (t.GetY() / divide), (int) (t.GetZ() / divide), t.GetVertexGlobalIds(0));
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
    auto key = std::make_pair(TMS_MakeGlobalVertexID(TMS_TrueParticles[i].GetRunID(), TMS_TrueParticles[i].GetVertexID()), TMS_TrueParticles[i].GetTrackId());
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

  EventCounter++;
}

TMS_Event::TMS_Event(TMS_Event &event, int slice) : TMS_Hits(event.GetHits(slice, true)), NonTMS_Hits(event.NonTMS_Hits),
      // Phase III: carry the truth side table over wholesale (simpler than filtering to just
      // this slice's hit IDs, and harmless -- lookups are still by HitId, which travels with
      // each TMS_Hit regardless of which subset ends up in this sliced event's TMS_Hits).
      HitIdCounter(event.HitIdCounter), TrueHitByHitId(event.TrueHitByHitId),
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
  RunNumber = event.RunNumber;
  SliceNumber = slice;
  SpillNumber = event.SpillNumber;
  
  
  nTrueTrajectories = -999;
  VertexIdOfMostEnergyInEvent = -9991;
  VertexGlobalIdOfMostEnergyInEvent = -9991;
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
    SetLeptonInfoUsingGlobalVertexID(GetVertexGlobalIdOfMostVisibleEnergy());
    long long primary_vertex_global_id = GetVertexGlobalIdOfMostVisibleEnergy();
    if (Reactions.find(primary_vertex_global_id) != Reactions.end())
      Reaction = Reactions[primary_vertex_global_id];
    else { Reaction = "NA"; std::cout<<"Warning: couldn't find reaction for primary vertex"<<std::endl; }
  }

  // Update the counts per slice
  ConnectTrueHitWithTrueParticle(true);
}

void TMS_Event::ApplyReconstructionEffects() {
  // First apply energy and timing models. Then merge hits. Then do a pedestal subtraction.
  // Sim-only steps (TMS_DetectorSimulation) and real-or-simulated steps (TMS_SignalProcessing)
  // are interleaved in this exact order deliberately: SimulateReadoutNoise() must run after
  // MergeCoincidentHits() so noise is drawn once per final merged channel-readout, not once
  // per raw sub-hit -- see TMS_DetectorSimulation::SimulateReadoutNoise().
  // Simulate an optical model
  TMS_DetectorSimulation::GetInstance().SimulateOpticalModel(*this, generator);
  // Noise hits can simulate hits
  TMS_DetectorSimulation::GetInstance().SimulateDarkCount(*this);
  // Simulate a timing model
  TMS_DetectorSimulation::GetInstance().SimulateTimingModel(*this, generator);
  // Simulate deadtime if needed
  TMS_DetectorSimulation::GetInstance().SimulateDeadtime(*this);
  // Merge hits that happened in the same scintillator strip and within the same readout time window
  TMS_SignalProcessing::GetInstance().MergeCoincidentHits(*this);
  // After merging hits, we have a single readout. This readout will have some electronic noise
  TMS_DetectorSimulation::GetInstance().SimulateReadoutNoise(*this, generator);
  // Simulate pedestal subtraction where any hit under Get_Sim_Readout_PedestalSubtractionThreshold is removed
  TMS_SignalProcessing::GetInstance().SimulatePedestalSubtraction(*this);
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
  for (const auto& it : Other_Event.TrueVisibleEnergyPerVertex) {
    TrueVisibleEnergyPerVertex[it.first] += it.second;
  }
  for (const auto& it : Other_Event.TrueVisibleEnergyPerParticle) {
    TrueVisibleEnergyPerParticle[it.first] += it.second;
  }
  for (const auto& it : Other_Event.Reactions) {
    auto inserted = Reactions.emplace(it.first, it.second);
    if (!inserted.second && inserted.first->second != it.second) {
      std::cout<<"Fatal: Conflicting reactions for global vertex id "<<it.first
               <<" while merging TMS events"<<std::endl;
      throw std::runtime_error("Fatal: conflicting reactions while merging events");
    }
  }
  for (const auto& it : Other_Event.info_about_vtx) {
    auto inserted = info_about_vtx.emplace(it.first, it.second);
    if (!inserted.second) {
      std::cout<<"Fatal: Duplicate global vertex id "<<it.first
               <<" while merging TMS events"<<std::endl;
      throw std::runtime_error("Fatal: duplicate global vertex id while merging events");
    }
  }
  // Reset this to recalculate on next call
  VertexIdOfMostEnergyInEvent = -9999;
  VertexGlobalIdOfMostEnergyInEvent = -9999;

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
  TrueNeutrinoPosition.SetX(vtx[0][0]*1000.);//change me to mm
  TrueNeutrinoPosition.SetY(vtx[0][1]*1000.);
  TrueNeutrinoPosition.SetZ(vtx[0][2]*1000.);
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
  // First find how much energy is in each variable
  for (auto& hit : TMS_Hits) {
    const TMS_TrueHit* true_hit = GetTrueHit(hit.GetHitId());
    if (true_hit == nullptr) continue; // No truth for this hit (e.g. real data)
    for (size_t i = 0; i < true_hit->GetNTrueParticles(); i++) {
      long long vertex_global_id = true_hit->GetVertexGlobalIds(i);
      // todo, true or reco energy?
      double energy = true_hit->GetEnergyShare(i);
      TrueVisibleEnergyPerVertex[vertex_global_id] += energy;
    }
  }
  
  // Now find the most energetic vertex
  double max = -1e9;
  long long max_vertex_global_id = -9992;
  double total_energy = 0;
  for (auto it : TrueVisibleEnergyPerVertex) {
    long long vertex = it.first;
    double energy = it.second;
    if (energy >= max) { max = energy; max_vertex_global_id = vertex; }
    total_energy += energy;
  }
  VertexGlobalIdOfMostEnergyInEvent = max_vertex_global_id;
  if (auto* vtx_info = GetVertexInfoByGlobalID(VertexGlobalIdOfMostEnergyInEvent)) {
    VertexIdOfMostEnergyInEvent = vtx_info->vtx_id;
  }
  else {
    VertexIdOfMostEnergyInEvent = -9992;
  }
  VisibleEnergyFromVertexInSlice = max;
  VisibleEnergyFromOtherVerticesInSlice = total_energy - max;
  
  if (TrueVisibleEnergyPerVertex.find(VertexGlobalIdOfMostEnergyInEvent) != TrueVisibleEnergyPerVertex.end())
    TotalVisibleEnergyFromVertex = TrueVisibleEnergyPerVertex[VertexGlobalIdOfMostEnergyInEvent];
  else {
    if (TrueVisibleEnergyPerVertex.size() > 0) {
      std::cout<<"Warning in GetVertexIdOfMostVisibleEnergy: TrueVisibleEnergyPerVertex.Find(VertexGlobalIdOfMostEnergyInEvent) == TrueVisibleEnergyPerVertex.end()"<<std::endl;
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

int TMS_Event::GetTrueParticleIndex(long long vertexglobalid, int trackid) {
  int out = -1;
  if (vertexglobalid >= 0 && trackid >= 0) {
    for (size_t i = 0; i < TMS_TrueParticles.size(); i++) {
      auto& tp = TMS_TrueParticles.at(i);
      if (TMS_MakeGlobalVertexID(tp.GetRunID(), tp.GetVertexID()) == vertexglobalid && tp.GetTrackId() == trackid) {
        out = i;
        break;
      }
    }
  }
  else {
    std::cout<<"GetTrueParticleIndex: Case of global vertex < 0 or trackid < 0. Global vertex id: "
             <<vertexglobalid<<", track id: "<<trackid
             <<", n TMS_TrueParticles: "<<TMS_TrueParticles.size()<<std::endl;
  }
  if (out < 0) std::cout<<"GetTrueParticleIndex: Case where out < 0. Global vertex id: "
                         <<vertexglobalid<<", track id: "<<trackid
                         <<", n TMS_TrueParticles: "<<TMS_TrueParticles.size()<<std::endl;
  return out;
}

int TMS_Event::GetPrimaryLeptonOfGlobalVertexID(long long vertexglobalid) {
  int lepton_index = -999;
  int current_index = 0;
  for (auto& particle : TMS_TruePrimaryParticles) {
    if (TMS_MakeGlobalVertexID(particle.GetRunID(), particle.GetVertexID()) == vertexglobalid) {
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

void TMS_Event::SetLeptonInfoUsingGlobalVertexID(long long vertexglobalid) {

  // And now fill lepton info
  auto particle_index = GetPrimaryLeptonOfGlobalVertexID(vertexglobalid);
  if (particle_index >= 0) {
    auto particle = TMS_TruePrimaryParticles[particle_index];
    int lepton_pdg = particle.GetPDG();
    auto lepton_position = particle.GetBirthPosition();
    auto lepton_momentum = particle.GetBirthMomentumAsLorentz();
    FillTrueLeptonInfo(lepton_pdg, lepton_position, lepton_momentum, particle.GetVertexID());
  }
  else {
    std::cout<<"Warning in SetLeptonInfoUsingGlobalVertexID: GetPrimaryLeptonOfGlobalVertexID didn't"
               "return a valid particle index for global vertex id "<<vertexglobalid<<std::endl;
    FillTrueLeptonInfo(-9999999, TLorentzVector(-9999999, -999999, -999999, -999999), 
      TLorentzVector(-9999999, -999999, -999999, -999999), -9999999);
  }
}

double TMS_Event::CalculateEnergyInLArOuterShell(double thickness, long long vertexglobalid) {
  double out = 0;
  // Lar doesn't have good timing info, so we want all non tms hits, not just in this slice
  for (const auto& hit : NonTMS_Hits) {
    if (vertexglobalid < 0 || hit.GetVertexGlobalIds(0) == vertexglobalid) {
      TVector3 position(hit.GetX(), hit.GetY(), hit.GetZ());
      if (TMS_Geom::GetInstance().IsInsideLAr(position) && !TMS_Geom::GetInstance().IsInsideLAr(position, thickness)) {
        out += hit.GetHadronicEnergy();
      }
    }
  }
  return out;
}

double TMS_Event::CalculateEnergyInLAr(long long vertexglobalid) {
  double out = 0;
  for (const auto& hit : NonTMS_Hits) {
    if (hit.GetVertexGlobalIds(0) < 0) std::cout<<"Warning: found true hit with < 0 VertexGlobalId"<<std::endl;
    if (vertexglobalid < 0 || hit.GetVertexGlobalIds(0) == vertexglobalid) {
      TVector3 position(hit.GetX(), hit.GetY(), hit.GetZ());
      if (TMS_Geom::GetInstance().IsInsideLAr(position))
        out += hit.GetHadronicEnergy();
    }
  }
  return out;
}


double TMS_Event::CalculateTotalNonTMSEnergy(long long vertexglobalid) {
  double out = 0;
  for (const auto& hit : NonTMS_Hits) {
    if (hit.GetVertexGlobalIds(0) < 0) std::cout<<"Warning: found true hit with < 0 VertexGlobalId"<<std::endl;
    if (vertexglobalid < 0 || hit.GetVertexGlobalIds(0) == vertexglobalid) out += hit.GetHadronicEnergy();
  }
  return out;
}

void TMS_Event::ConnectTrueHitWithTrueParticle(bool slice) {
  // Now count the number of true hits per particle
  std::map<std::pair<long long, int>, int> NHitsPerParticle;
  std::map<std::pair<long long, int>, double> EnergyPerParticle;
  for (auto& hit : TMS_Hits) {
    // Only count hits that are not ped subtracted
    if (!hit.GetPedSup()) {
      const TMS_TrueHit* true_hit = GetTrueHit(hit.GetHitId());
      if (true_hit == nullptr) continue; // No truth for this hit (e.g. real data)
      // Only add 1 hit for each key once, so track if we saw a key already
      std::map<std::pair<long long, int>, int> key_seen;
      for (size_t i = 0; i < true_hit->GetNTrueParticles(); i++) {
        auto key = std::make_pair(true_hit->GetVertexGlobalIds(i), true_hit->GetPrimaryIds(i));
        if (key_seen.find(key) == key_seen.end()) {
          NHitsPerParticle[key] += 1;
          key_seen[key] = 1;
        }
        EnergyPerParticle[key] += true_hit->GetEnergyShare(i);
      }
    }
  }
  for (size_t i = 0; i < TMS_TrueParticles.size(); i++) {
    int count = 0;
    double energy = 0;
    // If it's not in the map, don't create it
    auto key = std::make_pair(TMS_MakeGlobalVertexID(TMS_TrueParticles[i].GetRunID(), TMS_TrueParticles[i].GetVertexID()), TMS_TrueParticles[i].GetTrackId());
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
    const long long global_vertex_id = hit.GetVertexGlobalIds(i);
    if (info_about_vtx.find(global_vertex_id) != info_about_vtx.end()) {
      info_about_vtx[global_vertex_id].AddEnergyFromHit(hit, i);
    }
    else std::cout<<"This should not happen but I didn't find a vertex for global vertex id "
                  <<global_vertex_id<<std::endl;
  }
}

Vtx_Info* TMS_Event::GetVertexInfo(int run_id, int vertex_id) {
  Vtx_Info* out = NULL;
  const long long global_vertex_id = TMS_MakeGlobalVertexID(run_id, vertex_id);
  if (info_about_vtx.find(global_vertex_id) != info_about_vtx.end())
    out = &info_about_vtx.at(global_vertex_id);
  return out;
}

Vtx_Info* TMS_Event::GetVertexInfoByGlobalID(long long vertex_global_id) {
  Vtx_Info* out = NULL;
  if (info_about_vtx.find(vertex_global_id) != info_about_vtx.end())
    out = &info_about_vtx.at(vertex_global_id);
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
