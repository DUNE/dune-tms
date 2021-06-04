// For drawing the event
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TBox.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TROOT.h"

#include "TMS_Event.h"
#include "TMS_Reco.h"
#include "TMS_Manager.h"

class TMS_EventViewer {
  public:
    static TMS_EventViewer &GetViewer() {
      static TMS_EventViewer Instance;
      return Instance;
    }

    void Draw(TMS_Event &event);

    ~TMS_EventViewer() {
      if (nDraws > 0) {
        Canvas->Print(CanvasName+".pdf]");
        std::cout << "TMS_EventViewer drew " << nDraws << " events to " << CanvasName+".pdf" << std::endl;
      }
    }

  private:
    TMS_EventViewer();
    TMS_EventViewer(TMS_EventViewer const&) = delete;
    void operator=(TMS_EventViewer const&) = delete;

    TH2D *xz_view;
    TH2D *yz_view;
    TCanvas *Canvas;
    TString CanvasName;

    TBox *xz_box_FV;
    TBox *xz_box_Full;

    TBox *yz_box_FV;
    TBox *yz_box_Full;

    TBox *xz_dead_top;
    TBox *xz_dead_center;
    TBox *xz_dead_bottom;

    TLine *yz_Thin_Thick;
    TLine *xz_Thin_Thick;

    int nDraws;

    // Include track finding in draw
    bool DrawTrackFinding;
};
