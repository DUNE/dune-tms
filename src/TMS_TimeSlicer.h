#ifndef _TMS_TIMESLICER_H_SEEN_
#define _TMS_TIMESLICER_H_SEEN_
class TMS_Event;

class TMS_TimeSlicer {
  public:

    static TMS_TimeSlicer& GetSlicer() {
      static TMS_TimeSlicer Instance;
      return Instance;
    }

    int RunTimeSlicer(TMS_Event &event);
    int SimpleTimeSlicer(TMS_Event &event);
    // [Recon.Time] PerViewSlicing: slice each view on its own, then merge the
    // x- and y-view slices that overlap in time and z (see the .cpp).
    int PerViewTimeSlicer(TMS_Event &event);
    
};
#endif
