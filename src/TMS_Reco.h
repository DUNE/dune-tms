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

#include "TMS_Track.h"

#include "TMS_ChargeID.h"

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

      if      (Heuristic == HeuristicType::kManhattan) return deltax + deltay;
      else if (Heuristic == HeuristicType::kEuclidean) return sqrt(deltax * deltax + deltay * deltay);
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

// Maybe put this inside a separate namespace
class TMS_TrackFinder {

  public:

    static TMS_TrackFinder& GetFinder() {
      static TMS_TrackFinder Instance;
      return Instance;
    }

    void FindTracks(TMS_Event &event);
    const std::vector<TMS_Hit> & GetCandidatesU() { return CandidatesU; };
    const std::vector<TMS_Hit> & GetCandidatesV() { return CandidatesV; };
    const std::vector<TMS_Hit> & GetCandidatesX() { return CandidatesX; };
    const std::vector<std::vector<TMS_Hit> >& GetTotalCandidatesU() { return TotalCandidatesU; };
    const std::vector<std::vector<TMS_Hit> >& GetTotalCandidatesV() { return TotalCandidatesV; };
    const std::vector<std::vector<TMS_Hit> >& GetTotalCandidatesX() { return TotalCandidatesX; };


    const std::vector<TMS_Hit> &GetCleanedHits() { return CleanedHits; };

    TF1* GetHoughLineU() { return HoughLineU; };
    TF1* GetHoughLineV() { return HoughLineV; };
    TF1* GetHoughLineX() { return HoughLineX; };

    std::vector<std::pair<bool, TF1*> > GetHoughLinesU() { return HoughLinesU; };
    std::vector<std::pair<bool, TF1*> > GetHoughLinesV() { return HoughLinesV; };
    std::vector<std::pair<bool, TF1*> > GetHoughLinesX() { return HoughLinesX; };
    std::vector<std::pair<double,double> > GetHoughLinesU_Upstream() { return HoughLinesU_Upstream; };
    std::vector<std::pair<double,double> > GetHoughLinesU_Downstream() { return HoughLinesU_Downstream; };
    std::vector<std::pair<double,double> > GetHoughLinesV_Upstream() { return HoughLinesV_Upstream; };
    std::vector<std::pair<double,double> > GetHoughLinesV_Downstream() { return HoughLinesV_Downstream; };
    std::vector<std::pair<double,double> > GetHoughLinesX_Upstream() { return HoughLinesX_Upstream; };
    std::vector<std::pair<double,double> > GetHoughLinesX_Downstream() { return HoughLinesX_Downstream; };

    int **GetAccumulator() { return Accumulator; };

    std::vector<std::vector<TMS_Hit> > &GetClusterCandidatesU() { return ClusterCandidatesU; };
    std::vector<std::vector<TMS_Hit> > &GetClusterCandidatesV() { return ClusterCandidatesV; };
    std::vector<std::vector<TMS_Hit> > &GetClusterCandidatesX() { return ClusterCandidatesX; };
    std::vector<std::vector<TMS_Hit> > &GetHoughCandidatesU() { return HoughCandidatesU; };
    std::vector<std::vector<TMS_Hit> > &GetHoughCandidatesV() { return HoughCandidatesV; };
    std::vector<std::vector<TMS_Hit> > &GetHoughCandidatesX() { return HoughCandidatesX; };

    std::vector<TMS_Track> &GetHoughTracks3D() { return HoughTracks3D; };

    TH2D *AccumulatorToTH2D(bool zy);

    void SetZMaxHough(double z) { zMaxHough = z;};

    double CalculateTrackLength(const std::vector<TMS_Hit> &Hits);
    double CalculateTrackEnergy(const std::vector<TMS_Hit> &Hits);
    
    void CalculateTrackLengthU(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackEnergyU(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackLengthV(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackEnergyV(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackLengthX(const std::vector<std::vector<TMS_Hit> > &Hits);
    void CalculateTrackEnergyX(const std::vector<std::vector<TMS_Hit> > &Hits);

    double CalculateTrackLength3D(const TMS_Track &Hits);
    double CalculateTrackEnergy3D(const TMS_Track &Hits);
    double CalculateTrackKEByRange(const TMS_Track &Hits);

    double CalculateTrackLengthKalman(const TMS_Track &Hits);

    std::vector<std::vector<TMS_Hit> > FindClusters(const std::vector<TMS_Hit> &TMS_Hits);

    // Exclude hits Mask from set Orig
    void MaskHits(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask);
    void WalkDownStream(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask);
    void WalkUpStream(std::vector<TMS_Hit> &Orig, std::vector<TMS_Hit> &Mask);

    // Run a best first search
    void BestFirstSearch(const std::vector<TMS_Hit> &Hits, const char &hitgroup);

    //void HoughTransform(const std::vector<TMS_Hit> &Hits);
    std::vector<std::vector<TMS_Hit> > HoughTransform(const std::vector<TMS_Hit> &Hits, const char &hitgroup);
    std::vector<TMS_Hit> RunHough(const std::vector<TMS_Hit> &Hits, const char &hitgroup);

    std::vector<TMS_Hit> Extrapolation(const std::vector<TMS_Hit> &TrackHits, const std::vector<TMS_Hit> &Hits);
    std::vector<TMS_Track> TrackMatching3D();
    static bool SortByHitNumber(std::vector<TMS_Hit> &OneTrack, std::vector<TMS_Hit> &OtherTrack) { return (OneTrack.size() >= OtherTrack.size() ); };
    void FindPseudoXTrack();
    void CalculateRecoY(TMS_Hit &XHit, TMS_Hit &UHit, TMS_Hit &VHit); //XHit is here a placeholder for the hit that gets the RecoY set. This is then to be chosen as either UHit or VHit
    void CalculateRecoX(TMS_Hit &UHit, TMS_Hit &VHit, TMS_Hit &XHit);
    double CompareY(TMS_Hit &UHit, TMS_Hit &VHit, TMS_Hit &XHit);

    // Clean up the hits, removing duplicates and zero entries
    std::vector<TMS_Hit> CleanHits(const std::vector<TMS_Hit> &Hits);
    // Get hits projected onto xz or yz
    std::vector<TMS_Hit> ProjectHits(const std::vector<TMS_Hit> &Hits, TMS_Bar::BarType bartype = TMS_Bar::kXBar);
    std::vector<TMS_Hit> RunAstar(const std::vector<TMS_Hit> &Hits, bool ConnectAll = false);

    std::vector<TMS_Hit> UHitGroup;
    std::vector<TMS_Hit> VHitGroup;
    std::vector<TMS_Hit> XHitGroup;

    // Helper function to check if a hit is next to a gap
    bool NextToGap(double, double);

    void SpatialPrio(std::vector<TMS_Hit> &Hits);

    // Evaluate the track finding by using the event's true particles
    void EvaluateTrackFinding(TMS_Event &event);

    void ClearClass();

    std::vector<double> &GetTrackLengthU() { return TrackLengthU; };
    std::vector<double> &GetTrackLengthV() { return TrackLengthV; };
    std::vector<double> &GetTrackLengthX() { return TrackLengthX; };
    std::vector<double> &GetTrackEnergyU() { return TrackEnergyU; };
    std::vector<double> &GetTrackEnergyV() { return TrackEnergyV; };
    std::vector<double> &GetTrackEnergyX() { return TrackEnergyX; };

    void GetHoughLine(const std::vector<TMS_Hit> &TMS_Hits, double &slope, double &intercept);

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

    TMS_Kalman KalmanFilter;

    TMS_DBScan DBSCAN;
    TMS_ChargeID ChargeID;//ID_Track_Charge(const std::vector<TMS_Hit> &Hits);

    int FindBin(double Rho);
    // The candidates for each particle
    std::vector<TMS_Hit> CandidatesU;
    std::vector<TMS_Hit> CandidatesV;
    std::vector<TMS_Hit> CandidatesX;
    std::vector<TMS_Hit> RawHits;

    std::vector<TMS_Track> TotalTracks;
    std::vector<TMS_Hit> CleanedHits;
    std::vector<std::vector<TMS_Hit> > TotalCandidatesU;
    std::vector<std::vector<TMS_Hit> > TotalCandidatesV;
    std::vector<std::vector<TMS_Hit> > TotalCandidatesX;
    std::vector<std::pair<bool, TF1*> > HoughLinesU;
    std::vector<std::pair<bool, TF1*> > HoughLinesV;
    std::vector<std::pair<bool, TF1*> > HoughLinesX;
    std::vector<std::vector<TMS_Hit> > ClusterCandidatesU;
    std::vector<std::vector<TMS_Hit> > ClusterCandidatesV;
    std::vector<std::vector<TMS_Hit> > ClusterCandidatesX;
    std::vector<std::vector<TMS_Hit> > HoughCandidatesU;
    std::vector<std::vector<TMS_Hit> > HoughCandidatesV;
    std::vector<std::vector<TMS_Hit> > HoughCandidatesX;
    std::vector<std::vector<TMS_Hit> > SortedHoughCandidatesU;
    std::vector<std::vector<TMS_Hit> > SortedHoughCandidatesV;
    std::vector<std::vector<TMS_Hit> > SortedHoughCandidatesX;
    std::vector<std::pair<double,double>> HoughLinesU_Upstream;
    std::vector<std::pair<double,double>> HoughLinesU_Downstream;
    std::vector<std::pair<double,double>> HoughLinesV_Upstream;
    std::vector<std::pair<double,double>> HoughLinesV_Downstream;
    std::vector<std::pair<double,double>> HoughLinesX_Upstream;
    std::vector<std::pair<double,double>> HoughLinesX_Downstream;
    std::vector<TMS_Track> HoughTracks3D;

    char hitgroup;
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

    TF1 *HoughLineU;
    TF1 *HoughLineV;
    TF1 *HoughLineX;

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

    std::vector<double> TrackEnergyU;
    std::vector<double> TrackEnergyV;
    std::vector<double> TrackEnergyX;
    std::vector<double> TrackLengthU;
    std::vector<double> TrackLengthV;
    std::vector<double> TrackLengthX;

    // xvalue is x-axis, y value is y-axis
    void Accumulate(double xhit, double zhit); 
};

#endif
