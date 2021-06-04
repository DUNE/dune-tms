#ifndef __TMS_MANAGER_H__
#define __TMS_MANAGER_H__

// Just a global parameter manager
class TMS_Manager {

  public:
    static TMS_Manager& GetInstance() {
      static TMS_Manager Instance;
      return Instance;
    }

    void SetFileName(std::string file) { Filename = file; };
    std::string GetFileName() { return Filename; };

  private:
    TMS_Manager() {};
    TMS_Manager(TMS_Manager const &) = delete;
    void operator=(TMS_Manager const &) = delete;

    std::string Filename;
};

#endif
