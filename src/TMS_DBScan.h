#ifndef __TMS_DBSCAN_H__
#define __TMS_DBSCAN_H__

#include <cmath>
#include <vector>

// Enum for how a point has been classified
enum class PointClassifier { kUnclassified = -1, kNoise = 0 };

// Just a struct to keep the info
class TMS_DBScan_Point {
  public:
    TMS_DBScan_Point(double xvar, double yvar, double zvar) : x(xvar), y(yvar), z(zvar), ClusterID(int(PointClassifier::kUnclassified)) {};
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
      double xdist = x - PointTarget.x;
      double ydist = y - PointTarget.y;
      double zdist = z - PointTarget.z;
      return sqrt(xdist * xdist + ydist * ydist + zdist * zdist);
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
  {};
    TMS_DBScan() :
      //_Points(std::vector<TMS_DBScan_Point> Empty),
      _Epsilon(0),
      _MinPoints(0)
  {}

    // Some setters for the privates
    void SetPoints(std::vector<TMS_DBScan_Point> &Points) { _Points = Points; };
    void SetEpsilon(double eps) { 
      if (eps <= 0) {
        std::cerr << "Can't set epsilon in DBSCAN to negative or zero value!" << std::endl;
        throw;
      }
      _Epsilon = eps; 
    };
    void SetMinPoints(unsigned int minpts) { 
      if (minpts <= 0) {
        std::cerr << "Can't set minpts in DBSCAN to negative or zero value!" << std::endl;
        throw;
      }
      _MinPoints = minpts; 
    };

    // Return the points
    std::vector<TMS_DBScan_Point> &GetPoints() {return _Points;};

    void Print() {
      std::cout << "DBSCAN with " << _Points.size() << " points" << std::endl;
      std::cout << "   Found " << GetNClusters() << " clusters" << std::endl;
      std::cout << "         " << GetNNoise() << " noise" << std::endl;
      std::cout << "         " << GetNUnclassified() << " unclassified" << std::endl;
    }

    // Get the number of clusters
    int GetNClusters() {
      int HighestCluster = 0;
      for (auto &i : _Points) {
        if (i.ClusterID > HighestCluster) HighestCluster = i.ClusterID;
      }
      return HighestCluster;
    }

    // Get the clusters
    std::vector<std::vector<TMS_DBScan_Point> > GetClusters() {
      std::vector<std::vector<TMS_DBScan_Point> > vect;
      int nClusters = GetNClusters();
      vect.resize(nClusters);
      for (auto &i : _Points) {
        if (i.ClusterID > int(PointClassifier::kNoise)) {
          vect[i.ClusterID-1].push_back(std::move(i));
        }
      }
      return vect;
    }

    // Get all the unclassified noise points
    std::vector<TMS_DBScan_Point> GetNoise() {
      std::vector<TMS_DBScan_Point> vect;
      for (auto &i : _Points) {
        if (i.ClusterID == int(PointClassifier::kNoise)) {
          vect.push_back(std::move(i));
        }
      }
      return vect;
    }

    // Get the number of noise points
    int GetNNoise() {
      int NoiseHits = 0;
      for (auto &i : _Points) {
        if (i.ClusterID == int(PointClassifier::kNoise)) NoiseHits++;
      }
      return NoiseHits;
    }

    // Get the number of unclassified points
    int GetNUnclassified() {
      int Unclassified = 0;
      for (auto &i : _Points) {
        if (i.ClusterID == int(PointClassifier::kUnclassified)) Unclassified++;
      }
      return Unclassified;
    }

    // Run the scan
    bool RunDBScan() {
      if (_Points.empty()) {
#ifdef DEBUG
        std::cerr << "Can't run DBscan without having hits" << std::endl;
#endif
        return false;
      }
      // Start the cluster counter
      int ClusterID = 1; 
      // Start the loop over the points
      for (auto &i: _Points) {
        // Check if point has been classified already
        if (i.ClusterID != int(PointClassifier::kUnclassified)) continue;
          // This is the "inner loop"
          if (GrowCluster(i, ClusterID)) {
            ClusterID++;
          }
          // If this point has no neighbours it's noise
          else i.ClusterID = int(PointClassifier::kNoise);
      }
      // If we haven't found any clusters, the DBSCAN has failed
      if (ClusterID == 1) {
        return false;
      }
      return true;
    }

    bool GrowCluster(TMS_DBScan_Point &Point, int ClusterID) {
      // Find the neighbours
      std::vector<int> Neighbours = FindNeighbours(Point);
      // Check it has more neighbours than our minimum requirement
      if (Neighbours.size() < _MinPoints) {
        Point.ClusterID = int(PointClassifier::kNoise);
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
      Neighbours.erase(Neighbours.begin() + CorePointIndex);

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
          if (_Points[j].ClusterID == int(PointClassifier::kUnclassified) ||
              _Points[j].ClusterID == int(PointClassifier::kNoise)) {
            if (_Points[j].ClusterID == int(PointClassifier::kUnclassified)) {
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
    double _Epsilon;
    unsigned int _MinPoints; // The minimum number of points to make a cluster
};

#endif
