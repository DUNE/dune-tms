#ifndef __TMS_DBSCAN_H__
#define __TMS_DBSCAN_H__

#include <cmath>
#include <vector>

// Enum for how a point has been classified
enum PointClassifier { kUnclassified = -1, kNoise = 0, kCore, kBorder };

// Just a struct to keep the info
class TMS_DBScan_Point {
  public:
    TMS_DBScan_Point(double xvar, double yvar, double zvar) : x(xvar), y(yvar), z(zvar), ClusterID(kUnclassified) {};
    // Positions
    double x, y, z;
    // Which cluster does it belong to
    int ClusterID;
    // Utility function comparing positions
    bool IsSame(TMS_DBScan_Point &Other) {
      if (Other.x == x && Other.y == y && Other.z == z) return true;
      return false;
    }

    // Calculate the distance between two points
    double Distance(const TMS_DBScan_Point& PointTarget) {
      double xdist = x-PointTarget.x;
      double ydist = y-PointTarget.y;
      double zdist = z-PointTarget.z;
      return sqrt(xdist*xdist+ydist*ydist+zdist*zdist);
    }

    void Print() {
      std::cout << "TMS_DBScan_Point: " << std::endl;
      std::cout << "   {x,y,z} = {" << x << ", " << y << ", " << z << "}" << std::endl;
      std::cout << "   ClusterID = " << ClusterID << std::endl;
    }
};

// The main DBSCAN implementation
class TMS_DBScan {
  public:
    TMS_DBScan(unsigned int MinPoints, double Epsilon, std::vector<TMS_DBScan_Point> &Points) :
      _Points(Points),
      _Epsilon(Epsilon),
      _MinPoints(MinPoints)
  {
    std::cout << "Created DBSCAN instance with: " << std::endl;
    std::cout << "  MinPoints: " << _MinPoints << std::endl;
    std::cout << "  Points: " << _Points.size() << std::endl;
    std::cout << "  Epsilon: " << _Epsilon << std::endl;
  };

    std::vector<TMS_DBScan_Point> &GetPoints() {return _Points;};
    int GetNClusters() {
      int HighestCluster = 0;
      for (auto &i : _Points) {
        if (i.ClusterID > HighestCluster) HighestCluster = i.ClusterID;
      }
      return HighestCluster;
    }

    int GetNNoise() {
      int NoiseHits = 0;
      for (auto &i : _Points) {
        if (i.ClusterID == kNoise) NoiseHits++;
      }
      return NoiseHits;
    }

    int GetNUnclassified() {
      int Unclassified = 0;
      for (auto &i : _Points) {
        if (i.ClusterID == kUnclassified) Unclassified++;
      }
      return Unclassified;
    }

    void RunDBScan() {
      if (_Points.empty()) {
        std::cerr << "Can't run DBscan without having hits" << std::endl;
        return;
      }
      // Start the cluster counter
      int ClusterID = 1; 
      // Start the loop over the points
      for (auto &i: _Points) {
        // Check if point has been classified already
        if (i.ClusterID != kUnclassified) continue;
          // This is the "inner loop"
          if (GrowCluster(i, ClusterID)) {
            ClusterID++;
          }
          // If this point has no neighbours it's noise
          else i.ClusterID = kNoise;
      }
    }

    bool GrowCluster(TMS_DBScan_Point &Point, int ClusterID) {
      // Find the neighbours
      std::vector<int> Neighbours = FindNeighbours(Point);
      // Check it has more neighbours than our minimum requirement
      if (Neighbours.size() < _MinPoints) {
        Point.ClusterID = kNoise;
        return false;
      }

      int IndexIt = 0;
      int CorePointIndex = 0;
      // Loop over the neighbours
      for (auto &seed: Neighbours) {
        // Set all neighbours' cluster ID to be the same
        _Points[seed].ClusterID = ClusterID;
        // Find the Point we're actually scanning around
        if (_Points[seed].IsSame(Point)) {
          CorePointIndex = IndexIt;
        }
        IndexIt++;
      }

      // Remove the core point from the cluster
      Neighbours.erase(Neighbours.begin()+CorePointIndex);

      // Get the size of the neighbours of this point
      std::vector<int>::size_type n = Neighbours.size();

      // Now loop over the neighbours of this point (which is a neighbour of our main Point)
      for (std::vector<int>::size_type i = 0; i < n; ++i) {
        // Get the neighbours of this neighbour of this point
        std::vector<int> Neighbours_Neighbour = FindNeighbours(_Points[Neighbours[i]]);

        // If this neighbour has too few neighbours, move onto the next neighbour
        if (Neighbours_Neighbour.size() < _MinPoints) continue;

        // Go through to find the cluster IDs
        for (auto &j: Neighbours_Neighbour) {
          // If point is unclassified or noise, check it
          if (_Points[j].ClusterID == kUnclassified ||
              _Points[j].ClusterID == kNoise) {
            if (_Points[j].ClusterID == kUnclassified) {
              Neighbours.push_back(j);
              // Updated the n
              n = Neighbours.size();
            }
            _Points[j].ClusterID = ClusterID;
          }
        }
      }
      return true;
    }


    // Find the clusters to which this point belongs
    std::vector<int> FindNeighbours(TMS_DBScan_Point &Point) {
      // Loop over the clusters
      int NeighbourIndex = 0;
      // The clusters
      std::vector<int> NeighbourIndices;
      // Loop over the points
      for (auto &i: _Points) {
        if (Point.Distance(i) <= _Epsilon) NeighbourIndices.push_back(NeighbourIndex);
        NeighbourIndex++;
      }
      return NeighbourIndices;
    }


  private:
    std::vector<TMS_DBScan_Point> _Points; // The points
    unsigned int _MinPoints; // The minimum number of points to make a cluster
    double _Epsilon;
};

#endif
