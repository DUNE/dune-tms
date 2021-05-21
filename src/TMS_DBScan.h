#ifndef __TMS_DBSCAN_H__
#define __TMS_DBSCAN_H__

// Just a struct to keep the info
struct TMS_DBScan_Point {
  double x, y, z; // Positions
  int ClusterID; // Which cluster does it belong to
  // Don't want to check ClusterID
  bool IsSame(TMS_DBScan_Point &Other) {
    if (Other.x == x && Other.y == y && Other.z == z) return true;
    return false;
  }
}

// The main DBSCAN implementation
class TMS_DBScan {
  public:
    TMS_DBScan(unsigned int MinPoints, double Epsilon, std::vector<TMS_DBScan_Point> &Points) :
      _Points(Points),
      _Epsilon(Epsilon),
      _MinPoints(MinPoints)
  {
    }

    bool GrowCluster(TMS_DBScan_Point &Point, int ClusterID) {
      std::vector<int> ClusterSeeds = FindClusters(Point);
      if (ClusterSeeds.size() < _MinPoints) {
        Point.ClusterID = TMS_DBScan::kNoise;
        return false;
      }

      int IndexIt = 0;
      int CorePointIndex = 0;
      for (auto &seed: ClusterSeeds) {
        _Points[seed].ClusterID = ClusterID;
        // Check if point is the same
        if (_Points[seed].IsSame(Point)) {
          CorePointIndex = IndexIt;
        }
        IndexIt++;
      }
      // Remove the core point from the cluster
      ClusterSeeds.erase(ClusterSeeds.begin()+CorePointIndex);
    }

    // Calculate the distance between two points
    double Distance(const TMS_DBScan_Point& PointCore, const TMS_DBScan_Point& PointTarget) {
      double xdist = PointCore.x-PointTarget.x;
      double ydist = PointCore.y-PointTarget.y;
      double zdist = PointCore.z-PointTarget.z;
      return sqrt(xdist*xdist+ydist*ydist+zdist*zdist);
    }

    // Find the clusters to which this point belongs
    std::vector<int> FindClusters(TMS_DBScan_Point &Point) {
      // Loop over the clusters
      int ClusterIndex = 0;
      // The clusters
      std::vector<int> Clusters;
      // Loop over the points
      for (auto &i: _Points) {
        if (Distance(Point, *i) <= _Epsilon) Clusters.push_back(ClusterIndex);
        ClusterIndex++;
      }
      return Clusters;
    }

    // Enum for how a point has been classified
    enum PointClassifier { kUnclassified = 0, kCore, kBorder, kNoise };

  private:
    std::vector<TMS_DBScan_Point> _Points; // The points
    unsigned int _MinPoints; // The minimum number of points to make a cluster
    unsigned int _PointSize;
    double _Epsilon;

}

#endif
