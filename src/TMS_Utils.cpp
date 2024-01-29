#include "TMS_Utils.h"

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

  std::vector<std::pair<bool, TF1*> > HoughLinesOne = TMS_TrackFinder::GetFinder().GetHoughLinesOne();
  std::vector<std::pair<bool, TF1*> > HoughLinesOther = TMS_TrackFinder::GetFinder().GetHoughLinesOther();
  int nit = 0;
  for (auto it = HoughLinesOne.begin(); it != HoughLinesOne.end(); ++it, ++nit) {
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
    tracks[nit].TrackLength_gcm3 = TMS_TrackFinder::GetFinder().GetTrackLengthOne()[nit];
    tracks[nit].TrackEnergy = TMS_TrackFinder::GetFinder().GetTrackEnergyOne()[nit];
  }
  keeper_nit = nit;
  nit = 0;
  for (auto it = HoughLinesOther.begin(); it != HoughLinesOther.end(); ++it, ++nit)  {
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
    tracks[nit+keeper_nit].TrackLength_gcm3 = TMS_TrackFinder::GetFinder().GetTrackLengthOther()[nit];
    tracks[nit+keeper_nit].TrackEnergy = TMS_TrackFinder::GetFinder().GetTrackEnergyOther()[nit];
  }

  return caf;
}
#endif
