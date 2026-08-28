#ifndef _TMS_SPACEPOINT_H_SEEN_
#define _TMS_SPACEPOINT_H_SEEN_

class TMS_SpacePoint {
  public:
    TMS_SpacePoint() :
      x(0), y(0), z(0),
      x_hit_index(-1), y_hit_index(-1),
      time(0) {}

    TMS_SpacePoint(double x_val, double y_val, double z_val,
                   int x_idx, int y_idx, double time_val) :
      x(x_val), y(y_val), z(z_val),
      x_hit_index(x_idx), y_hit_index(y_idx),
      time(time_val) {}

    // Getters
    double GetX() const { return x; }
    double GetY() const { return y; }
    double GetZ() const { return z; }
    int GetXHitIndex() const { return x_hit_index; }
    int GetYHitIndex() const { return y_hit_index; }
    double GetTime() const { return time; }

    // Setters
    void SetPosition(double x_val, double y_val, double z_val) {
      x = x_val;
      y = y_val;
      z = z_val;
    }
    void SetHitIndices(int x_idx, int y_idx) {
      x_hit_index = x_idx;
      y_hit_index = y_idx;
    }
    void SetTime(double time_val) { time = time_val; }

  private:
    double x;
    double y;
    double z;
    int x_hit_index;
    int y_hit_index;
    double time;
};

#endif
