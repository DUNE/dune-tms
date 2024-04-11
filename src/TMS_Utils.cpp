#include "TMS_Utils.h"

#include "TMS_Reco.h"



//caf::SRTMS TMS_CAF_converter::ConvertEvent(TMS_Event &event) {
#ifdef DUNEANAOBJ_ENABLED
caf::SRTMS ConvertEvent() {

  std::vector<std::vector<TMS_Hit> > LineCandidates = TMS_TrackFinder::GetFinder().GetLineCandidates();

  // CAF object for TMS
  caf::SRTMS caf;
  caf.ntracks = LineCandidates.size();
  // Loop over tracks and convert them to caf::SRTrack format
  std::vector<caf::SRTrack> tracks = caf.tracks;
  tracks.reserve(caf.ntracks);

  // Build up number of total hits
  int TotalHits = TMS_TrackFinder::GetFinder().GetCleanedHits().size();

  // Fill the number of tracks from the track finder
  for (auto it = LineCandidates.begin(); it != LineCandidates.end(); ++it) {
    // Get the track
    auto track = *it;

    // caf tracks are pretty simple right now; 
    // only start point, end point, start, direction and end direction
    TMS_Hit firsthit = track.front();
    TMS_Hit lasthit = track.back();
    // Convert the TMS_Hit hit information to SRVector3D
    caf::SRVector3D firsthit_conv(firsthit.GetNotZ(), 200, firsthit.GetZ()); // Have no y info in the TMS
    caf::SRVector3D lasthit_conv(lasthit.GetNotZ(), 200, lasthit.GetZ()); // Have no y info in the TMS

    caf::SRTrack srtrack;
    srtrack.start = firsthit_conv;
    srtrack.end = lasthit_conv;
    srtrack.TrackQuality = double((*it).size())/TotalHits;
    caf.tracks.push_back(srtrack);
  }

  std::vector<std::pair<bool, TF1*> > HoughLinesU = TMS_TrackFinder::GetFinder().GetHoughLinesU();
  std::vector<std::pair<bool, TF1*> > HoughLinesV = TMS_TrackFinder::GetFinder().GetHoughLinesV();
  std::vector<std::pair<bool, TF1*> > HoughLinesX = TMS_TrackFinder::GetFinder().GetHoughLinesX();
  int nit = 0;
  for (auto it = HoughLinesU.begin(); it != HoughLinesU.end(); ++it, ++nit) {
    //double intercept = (*it).second->GetParameter(0);
    //double slope = (*it).second->GetParameter(1);

    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double zlow = TMS_Const::TMS_Thin_Start;
    double zhi = TMS_Const::TMS_Thick_Start;
    double xlow = (*it).second->Eval(zlow);
    double xhi = (*it).second->Eval(zhi);

    double zlen = zhi-zlow;
    double xlen = xhi-xlow;
    double len = sqrt(xlen*xlen+zlen*zlen);
    zlen = zlen/len;
    xlen = xlen/len;

    caf::SRVector3D dir(xlen, 0, zlen); // Make the converted direction unit vector
    tracks[nit].dir = dir;

    // Now do the track quality, track energy and track length
    tracks[nit].TrackLength_gcm3 = TMS_TrackFinder::GetFinder().GetTrackLengthU()[nit];
    tracks[nit].TrackEnergy = TMS_TrackFinder::GetFinder().GetTrackEnergyU()[nit];
  }
  keeper_nit = nit;
  nit = 0;
  for (auto it = HoughLinesV.begin(); it != HoughLinesV.end(); ++it, ++nit)  {
    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double zlow = TMS_Const::TMS_Thin_Start;
    double zhi = TMS_Const::TMS_Thick_Start;
    double xlow = (*it).second->Eval(zlow);
    double xhi = (*it).second->Eval(zhi);

    double zlen = zhi-zlow;
    double xlen = xhi-xlow;
    double len = sqrt(xlen*xlen+zlen*zlen);
    zlen = zlen/len;
    xlen = xlen/len;

    caf::SRVector3D dir(xlen, 0, zlen); // Make the converted direction unit vector
    tracks[nit+keeper_nit].dir = dir;

    // Now do the track quality, track energy and track length
    tracks[nit+keeper_nit].TrackLength_gcm3 = TMS_TrackFinder::GetFinder().GetTrackLengthV()[nit];
    tracks[nit+keeper_nit].TrackEnergy = TMS_TrackFinder::GetFinder().GetTrackEnergyV()[nit];
  }

  keeper_nit += nit;
  nit = 0;
  for (auto it = HoughLinesX.begin(); it != HoughLinesX.end(); ++it, ++nit)  {
    // Calculate the z and x vectors by evaling the TF1 in thin and thick target
    double zlow = TMS_Const::TMS_Thin_Start;
    double zhi = TMS_Const::TMS_Thick_Start;
    double xlow = (*it).second->Eval(zlow);
    double xhi = (*it).second->Eval(zhi);

    double zlen = zhi-zlow;
    double xlen = xhi-xlow;
    double len = sqrt(xlen*xlen+zlen*zlen);
    zlen = zlen/len;
    xlen = xlen/len;

    caf::SRVector3D dir(xlen, 0, zlen); // Make the converted direction unit vector
    tracks[nit+keeper_nit].dir = dir;

    // Now do the track quality, track energy and track length
    tracks[nit+keeper_nit].TrackLength_gcm3 = TMS_TrackFinder::GetFinder().GetTrackLengthX()[nit];
    tracks[nit+keeper_nit].TrackEnergy = TMS_TrackFinder::GetFinder().GetTrackEnergyX()[nit];
  }
 

  return caf;
}
#endif


namespace TMS_Utils {
  TMS_Utils::ParticleInfo GetPrimaryIdsByEnergy(const std::vector<TMS_Hit>& hits) { 
      std::unordered_map<int, double> totalMap;

      // Iterate through the list of hits
      int total_n_true_particles = 0;
      for (const auto& hit : hits) {
        auto true_hit = hit.GetTrueHit();
        total_n_true_particles += true_hit.GetNTrueParticles();
        for (size_t i = 0; i < true_hit.GetNTrueParticles(); i++) {
          int pid = true_hit.GetPrimaryIds(i);
          double energy = true_hit.GetEnergyShare(i);
          // Add the utility to the corresponding pid in the map
          totalMap[pid] += energy;
        }
      }

      auto out = TMS_Utils::GetSumAndHighest(totalMap);
      // todo, remove when dark noise is added
      if (total_n_true_particles == 0 && hits.size() > 0) {
        std::cout<<"Warning: Did not find true particles in GetPrimaryIdsByEnergy.\n";
        std::cout<<"This could only happen with dark noise but that's not a thing yet\n";
        std::cout<<"Starting with n hits: "<<hits.size()<<"\n";
        std::cout<<"that contain this many n true particles: "<<total_n_true_particles<<"\n";
        std::cout<<"Found total energy: "<<out.total_energy<<"\n";
        std::cout<<"Found n indices: "<<out.indices.size()<<"\n";
        std::cout<<"Found total energy: "<<out.energies.size()<<std::endl;
      }
      return out;
  }

  struct ParticlePair {
    int index;
    double energy;
    
    ParticlePair(int i, double e) : index(i), energy(e) {}
  };

  static bool CompareByEnergy(const ParticlePair& a, const ParticlePair& b) {
      return a.energy > b.energy;
  }


  TMS_Utils::ParticleInfo GetSumAndHighest(const std::unordered_map<int, double>& map) {
      // First make a vector and sort it
      std::vector<ParticlePair> pairs;
      for (const auto& pair : map) {
          pairs.push_back(ParticlePair(pair.first, pair.second));
      }
      std::sort(pairs.begin(), pairs.end(), CompareByEnergy);

      TMS_Utils::ParticleInfo out;
      for (const auto& pair : pairs) {
        out.total_energy += pair.energy;
        out.indices.push_back(pair.index);
        out.energies.push_back(pair.energy);
      }
      
      return out;
  }

}




