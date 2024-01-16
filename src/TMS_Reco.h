#ifndef _TMS_RECO_H_SEEN_
#define _TMS_RECO_H_SEEN_

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <queue>
#include <list>
#include <unordered_map>

// For Hough line visualisation
#include "TF1.h"
#include "TH2D.h"

// TMS includes
#include "TMS_Manager.h"

#include "TMS_Hit.h"
#include "TMS_Event.h"
#include "TMS_Constants.h"
#include "TMS_DBScan.h"
#include "TMS_Utils.h"

// Hand over to the Kalman reconstruction once we find tracks
#include "TMS_Kalman.h"

#define __LARGE_COST__ 999999999

enum class TrackMethod {
  kHough,
  kAStar,
  kDBSCAN,
  kUnknown
};

enum class HeuristicType { 
  kManhattan, 
  kEuclidean, 
  kDetectorZ,
  kDetectorNotZ,
  kUnknown
};

// Utility class struct to store the node for track finding using A* or Best-First
class aNode {
  public:

    //aNode(double xval, double yval, double ywval): 
    aNode(double xval, double yval) :
      x(xval), y(yval),
      HeuristicCost(__LARGE_COST__), NodeID(-999),
      //Heuristic(kManhattan) { // what calculator
      Heuristic(HeuristicType::kEuclidean) { // what calculator
    };

    aNode(double xval, double yval, int ID): aNode(xval, yval) {
      NodeID = ID;
    };

    bool operator==(aNode const &other) {
      return (x == other.x && y == other.y);
    }

    double CalculateGroundCost(aNode const &other) {
      double deltax = std::abs(x-other.x);
      double deltay = std::abs(y-other.y);

      // Allow for all connections that exceed 1 cell, but if they're not great make sure they're heavily penalised
      if (deltay > 1 || deltax > 1) {
        return (deltax*deltax+deltay*deltay)*__LARGE_COST__;
      }

      double GroundCost = 0;
      // First add up the individual distance
      //GroundCost += (deltax+2*deltay)*10;
      GroundCost += (deltax*deltax + deltay*deltay)*(deltax*deltax + deltay*deltay);

      //if      (deltax+deltay == 2) GroundCost += 20;
      //if      (deltax+deltay == 2) GroundCost += 5;

      // Need to penalise diagonal connections to avoid them being preferred over non-diagonal
//      if (deltay == 1 && deltax == 1) GroundCost += 100;
      //if (deltax == 2 || deltay == 2) {
        //GroundCost += 10;
        //if (deltax == 2 && deltay == 2) {
          //GroundCost += 10;
        //}
      //}

      //else if (deltax+deltay == 3) GroundCost += 30;
      //else if (deltax+deltay == 4) GroundCost += 40;
      //else if (deltax+deltay == 5) GroundCost += 50;

      // Penalise double jumps
      //if (deltax > 1) GroundCost *= 10;

      /*
      if      (deltax == 1) GroundCost += deltax*10;
      else if (deltax == 2) GroundCost += deltax*100;

      if      (deltay == 1) GroundCost += deltay*10;
      else if (deltay == 2) GroundCost += deltay*100;
      else if (deltay == 3) GroundCost += deltay*200;

      if      (deltax == 1 && deltay == 1) GroundCost += 10;
      else if (deltax == 1 && deltay == 2) GroundCost += 100;
      else if (deltax == 1 && deltay == 3) GroundCost += 1000;

      if      (deltax == 2 && deltay == 1) GroundCost += 10;
      else if (deltax == 2 && deltay == 2) GroundCost += 100;
      else if (deltax == 2 && deltay == 3) GroundCost += 1000;
      */

      return GroundCost;
    }

    double Calculate(aNode const &other) {
      // x is the plane number, y is in mm
      // jumping one plane incurs 10 ground, so reflect that here; jumping 2 planes (i.e. adjacent) should be 10 ground, jumping 4 planes (i.e. next to adjacent) is double that
      //double deltax = (x-other.x)*5;
      double deltax = std::abs((x-other.x)*1);  //10->1
      // Moving 1 'bar' up is 10 ground cost, so reflect that here too (the bar width is here assumed to be 1cm)
      double deltay = std::abs((y-other.y)*1);  //10->1

      if      (Heuristic == HeuristicType::kManhattan) return deltax+deltay;
      else if (Heuristic == HeuristicType::kEuclidean) return sqrt(deltax*deltax+deltay*deltay);
      else if (Heuristic == HeuristicType::kDetectorZ) return deltax;
      else if (Heuristic == HeuristicType::kDetectorNotZ) return deltay;

      return __LARGE_COST__;
    }

    void SetHeuristicCost(aNode const &other) {
      HeuristicCost = Calculate(other);
    }

    void SetHeuristicCost(double val) {
      HeuristicCost = val;
    }

    void SetHeuristic(HeuristicType a) {
      Heuristic = a;
    }

    // Position
    double x;
    double y;
    // Costs
    double HeuristicCost;
    // Neighbours
    std::unordered_map<aNode*, double> Neighbours;
    // self ID
    int NodeID;
    // What Heuristic do we use
    HeuristicType Heuristic;

    void Print() const {
      std::cout << "NodeID: " << NodeID << std::endl;
      std::cout << "x, y = " << x << ", " << y << std::endl;
      std::cout << "Heuristic: " << HeuristicCost << std::endl; //" Ground: " << GroundCost << std::endl;
      std::cout << "Number of neighbours: " << Neighbours.size() << std::endl;
    }
};

inline bool operator<(aNode const &a, aNode const &b) {
  return a.HeuristicCost < b.HeuristicCost;
}

class TMS_TimeSlicer {
  public:

    static TMS_TimeSlicer& GetSlicer() {
      static TMS_TimeSlicer Instance;
      return Instance;
    }

    int RunTimeSlicer(TMS_Event &event);
    int SimpleTimeSlicer(TMS_Event &event);
    
};

// Maybe put this inside a separate namespace
class TMS_TrackFinder {

  public:

    static TMS_TrackFinder& GetFinder() {
      static TMS_TrackFinder Instance;
      return Instance;
    }

    void FindTracks(TMS_Event &event);
    const std::vector<TMS_Hit> & GetCandidatesOne() { return CandidatesOne; };
    const std::vector<TMS_Hit> & GetCandidatesOther() { return CandidatesOther; };
    const std::vector<std::vector<TMS_Hit> >& GetTotalCandidatesOne() { return TotalCandidatesOne; };
    const std::vector<std::vector<TMS_Hit> >& GetTotalCandidatesOther() { return TotalCandidatesOther; };

    const std::vector<TMS_Hit> &GetCleanedHits() { return CleanedHits; };

    TF1* GetHoughLineOne() { return HoughLineOne; };
    TF1* GetHoughLineOther() { return HoughLineOther; };

    std::vector<std::pair<bool, TF1*> > GetHoughLinesOne() { return HoughLinesOne; };
    std::vector<std::pair<bool, TF1*> > GetHoughLinesOther() { return HoughLinesOther; };
    std::vector<std::pair<double,double> > GetHoughLinesOne_Upstream() { return HoughLinesOne_Upstream; };
    std::vector<std::pair<double,double> > GetHoughLinesOne_Downstream() { return HoughLinesOne_Downstream; };
    std::vector<std::pair<double,double> > GetHoughLinesOther_Upstream() { return HoughLinesOther_Upstream; };
    std::vector<std::pair<double,double> > GetHoughLinesOther_Downstream() { return HoughLinesOther_Downstream; };

    int **GetAccumulator() { return Accumulator; };

    std::vector<std::vector<TMS_Hit> > &GetClusterCandidatesOne() { return ClusterCandidatesOne; };
    std::vector<std::vector<TMS_Hit> > &GetClusterCandidatesOther() { return ClusterCandidatesOther; };
    std::vector<std::vector<TMS_Hit> > &GetHoughCandidatesOne() { return HoughCandidatesOne; };
    std::vector<std::vector<TMS_Hit> > &GetHoughCandidatesOther() { return HoughCandidatesOther; };

    std::vector<std::vector<TMS_Hit> > &GetHoughTrack3D() { return HoughTrack3D; };

    TH2D *AccumulatorToTH2D(bool zy);

    void SetZMaxHough(double z) { zMaxHough = z;};

    void CalculateTrackLengthOne();
    void CalculateTrackLengthOther();
    void CalculateTrackEnergyOne();
    void CalculateTrackEnergyOther();
    
    void CalculateTrackLengthOne(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackEnergyOne(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackLengthOther(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackEnergyOther(const std::vector<std::vector<TMS_Hit> > &Hits);

    void CalculateTrackLength3D();
    void CalculateTrackEnergy3D();

    std::vector<std::vector<TMS_Hit> > FindClusters(const std::vector<TMS_Hit> &TMS_Hits);

    // Exclude hits Mask from set Orig
    void MaskHits(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask);
    void WalkDownStream(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask);
    void WalkUpStream(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask);

    // Run a best first search
    void BestFirstSearch(const std::vector<TMS_Hit> &Hits, const int &hitgroup);

    //void HoughTransform(const std::vector<TMS_Hit> &Hits);
    std::vector<std::vector<TMS_Hit> > HoughTransform(const std::vector<TMS_Hit> &Hits, const int &hitgroup);
    std::vector<TMS_Hit> RunHough(const std::vector<TMS_Hit> &Hits, const int &hitgroup);

    std::vector<TMS_Hit> Extrapolation(const std::vector<TMS_Hit> &TrackHits, const std::vector<TMS_Hit> &Hits);
    std::vector<std::vector<TMS_Hit> > TrackMatching3D();

    // Clean up the hits, removing duplicates and zero entries
    std::vector<TMS_Hit> CleanHits(const std::vector<TMS_Hit> &Hits);
    // Get hits projected onto xz or yz
    std::vector<TMS_Hit> ProjectHits(const std::vector<TMS_Hit> &Hits, TMS_Bar::BarType bartype = TMS_Bar::kXBar);
    std::vector<TMS_Hit> RunAstar(const std::vector<TMS_Hit> &Hits, bool ConnectAll = false);

    std::vector<TMS_Hit> OneHitGroup;
    std::vector<TMS_Hit> OtherHitGroup;

    // Helper function to check if a hit is next to a gap
    bool NextToGap(double, double);

    void SpatialPrio(std::vector<TMS_Hit> &Hits);

    // Evaluate the track finding by using the event's true particles
    void EvaluateTrackFinding(TMS_Event &event);

    void ClearClass();

    std::vector<double> &GetTrackLengthOne() { return TrackLengthOne; };
    std::vector<double> &GetTrackLengthOther() { return TrackLengthOther; };
    std::vector<double> &GetTrackEnergyOne() { return TrackEnergyOne; };
    std::vector<double> &GetTrackEnergyOther() { return TrackEnergyOther; };

    std::vector<double> &GetTrackLength3D() { return TrackLength3D; };
    std::vector<double> &GetTrackEnergy3D() { return TrackEnergy3D; };

    void GetHoughLine(const std::vector<TMS_Hit> &TMS_Hits, double &slope, double &intercept) {
      // Reset the accumulator
      for (int i = 0; i < nSlope; ++i) {
        for (int j = 0; j < nIntercept; ++j) {
          Accumulator[i][j] = 0;
        }
      }
      
      // First run a simple Hough Transform
      for (std::vector<TMS_Hit>::const_iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
        TMS_Hit hit = (*it);
        double xhit = hit.GetNotZ();
        double zhit = hit.GetZ();
        
        // If z position is above region of interest, ignore hit
        //if (IsXZ && zhit > zMaxHough) continue;
        if (zhit > zMaxHough) continue;
        
        Accumulate(xhit, zhit);
        
      }
      
      // Find the maximum of the accumulator and which m,c bin the maximum occurs in
      double max_zy = 0;
      int max_zy_slope_bin = 0;
      int max_zy_inter_bin = 0;
      for (int i = 0; i < nSlope; ++i) {
        for (int j = 0; j < nIntercept; ++j) {
          if (Accumulator[i][j] > max_zy) {
            max_zy = Accumulator[i][j];
            max_zy_slope_bin = i;
            max_zy_inter_bin = j;
          }
        }
      }
      intercept = InterceptMin+max_zy_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
      slope = SlopeMin+max_zy_slope_bin*(SlopeMax-SlopeMin)/nSlope;
    }


    TH1D* GetEfficiencyHist() { return Efficiency; };
    TH1D* GetTotalHist() { return Total; };
    TH1D* GetEfficiency() { 
      TH1D *eff = (TH1D*)Efficiency->Clone("Efficiency_ratio");
      eff->Divide(Total);
      return eff;
    }

  private:
    TMS_TrackFinder();
    TMS_TrackFinder(TMS_TrackFinder const &) = delete;
    void operator=(TMS_TrackFinder const &) = delete;
    ~TMS_TrackFinder() {};

    TMS_Kalman KalmanFitter;
    TMS_DBScan DBSCAN;

    int FindBin(double Rho);
    // The candidates for each particle
    std::vector<TMS_Hit> CandidatesOne;
    std::vector<TMS_Hit> CandidatesOther;
    std::vector<TMS_Hit> RawHits;

    std::vector<TMS_Hit> CleanedHits;
    std::vector<std::vector<TMS_Hit> > TotalCandidatesOne;
    std::vector<std::vector<TMS_Hit> > TotalCandidatesOther;
    std::vector<std::pair<bool, TF1*> > HoughLinesOne;
    std::vector<std::pair<bool, TF1*> > HoughLinesOther;
    std::vector<std::vector<TMS_Hit> > ClusterCandidatesOne;
    std::vector<std::vector<TMS_Hit> > ClusterCandidatesOther;
    std::vector<std::vector<TMS_Hit> > HoughCandidatesOne;
    std::vector<std::vector<TMS_Hit> > HoughCandidatesOther;
    std::vector<std::pair<double,double>> HoughLinesOne_Upstream;
    std::vector<std::pair<double,double>> HoughLinesOne_Downstream;
    std::vector<std::pair<double,double>> HoughLinesOther_Upstream;
    std::vector<std::pair<double,double>> HoughLinesOther_Downstream;
    std::vector<std::vector<TMS_Hit> > HoughTrack3D;

    int hitgroup;

    int nIntercept;
    int nSlope;
    double InterceptMin;
    double InterceptMax;
    double SlopeMin;
    double SlopeMax;
    double InterceptWidth;
    double SlopeWidth;

    // Maximum that we run Hough in z
    double zMinHough;
    double zMaxHough;
    // Maximum number of Hough transforms
    int nMaxHough;
    // How many continous hits in the first layers to use Hough in thin layer only
    int nThinCont;
    // Number of hits required to actually run the Hough transform
    double nHits_Tol;

    int **Accumulator;

    TF1 *HoughLineOne;
    TF1 *HoughLineOther;

    unsigned int nMinHits;
    unsigned int nMaxMerges;
    double MinDistHough;

    bool IsGreedy;
    // Which planes are next to the gaps (i.e. may cause discontinuities)?
    std::vector<int> PlanesNearGap;

    TH1D *Efficiency;
    TH1D *Total;

    HeuristicType kHeuristic;

    bool UseClustering;
    TrackMethod kTrackMethod;

    std::vector<double> TrackEnergyOne;
    std::vector<double> TrackEnergyOther;
    std::vector<double> TrackLengthOne;
    std::vector<double> TrackLengthOther;

    std::vector<double> TrackEnergy3D;
    std::vector<double> TrackLength3D;

    // xvalue is x-axis, y value is y-axis
    void Accumulate(double xhit, double zhit) {
      
      // Could probably multi-thread this operation
      // Now do the Hough
      for (int i = 0; i < nSlope; ++i) {
        double m = SlopeMin+i*SlopeWidth;
        if (m > SlopeMax) m = SlopeMax;
        
        // Now calculate rho
        double c = xhit-m*zhit;
        if (c > InterceptMax) c = InterceptMax;
        
        // Find which rho bin this corresponds to
        int c_bin = FindBin(c);
        
        /*
           if (i > nSlope || c_bin > nIntercept) {
           std::cout << "c: " << c << std::endl;
           std::cout << "m: " << m << std::endl;
           std::cout << "i: " <<  i << std::endl;
           std::cout << "cbin: " << c_bin << std::endl;
           }
           */

        // Fill the accumulator
        Accumulator[i][c_bin]++;
      }
    }
};

#endif
