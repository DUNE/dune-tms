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
    
};
#endif
