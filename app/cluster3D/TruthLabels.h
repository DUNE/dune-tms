#ifndef _TRUTHLABELS_H_SEEN_
#define _TRUTHLABELS_H_SEEN_

#include <cstddef>
#include <functional>

// Truth bookkeeping shared by the Cluster3D validation tools.

// A true particle: vertex global id and (collapsed) track id.
struct TrueLabel {
  long long vgid = -1;
  int trackid = -999;
  bool Valid() const { return vgid >= 0; }
  bool operator==(const TrueLabel &o) const { return vgid == o.vgid && trackid == o.trackid; }
};

struct LabelHash {
  std::size_t operator()(const TrueLabel &l) const {
    return std::hash<long long>()(l.vgid) ^ (std::hash<int>()(l.trackid) << 1);
  }
};

// Per-hit truth: the hit's two largest true contributors (primaries, with
// their own secondaries folded in) and their energy fractions, from the
// SpacePointHitTrue* branches. A hit two particles cross is genuinely shared,
// so metrics credit each its share.
struct HitTruth {
  TrueLabel first, second;
  double first_frac = 0.0, second_frac = 0.0;
  double Share(const TrueLabel &who) const {
    return (first.Valid() && first == who ? first_frac : 0.0) +
           (second.Valid() && second == who ? second_frac : 0.0);
  }
};

// Owner of a space point from its two hits' truth: the particle with the
// largest mean energy share over the two hits, if that mean exceeds 0.5 --
// the point has to be more than half one particle's. A genuine point scores
// 1; a ghost pairing one particle's hit with another's scores 0.5 and has no
// owner (2026-09-25: case G's fit rode a ghost track -- one muon's x hits
// with another muon's y hits -- which single-sided labels credited 100%).
inline TrueLabel PointOwner(const HitTruth &a, const HitTruth &b) {
  TrueLabel best;
  double bestScore = 0.5;
  for (const TrueLabel &who : {a.first, a.second, b.first, b.second}) {
    if (!who.Valid()) continue;
    const double score = 0.5 * (a.Share(who) + b.Share(who));
    if (score > bestScore) {
      bestScore = score;
      best = who;
    }
  }
  return best;
}

#endif
