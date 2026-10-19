#ifndef _TMS_BADCHANNELS_H_SEEN_
#define _TMS_BADCHANNELS_H_SEEN_

#include <map>
#include <string>

#include "TMS_Hit.h"

// Electronics channels (a bar, or a bar half for the split X-bars: TMS_ChannelId) that are known
// not to give a usable readout. The decoder drops their readouts before any TMS_Hit is made, so
// they are applied the same way to simulation and, later, to data.
//
// The list is read from the file named by [Decoder] BadChannelFile in the main TOML. Each entry:
//   [[BadChannel]]
//   Plane = 12, Bar = 40, View = 3   # TMS_ChannelId; bar numbers are unique only within a plane
//   Status = "dead"                  # "missing" (no such bar) or "dead" (bar does not report)
//   Reason = "free text"             # optional
// A status other than those two is an error, so a misspelled status cannot silently pass as good.
// The shipped list is empty. A database with validity ranges can replace the file behind this class
// without changing how the decoder asks (IsBad).
class TMS_BadChannels {
  public:
    enum class Status { Good, Missing, Dead };

    struct Entry {
      Status State;
      std::string Reason;
    };

    // Reads the file; "" means no bad channels. A relative path is taken from $TMS_DIR.
    // Throws if the file cannot be read or an entry is malformed.
    explicit TMS_BadChannels(const std::string& filename = "");

    // The list named in the main TOML, read on first use
    static const TMS_BadChannels& GetInstance();

    Status GetStatus(const TMS_ChannelId& channel) const;
    // True for every state whose readouts are dropped (Missing and Dead for now)
    bool IsBad(const TMS_ChannelId& channel) const { return GetStatus(channel) != Status::Good; }
    size_t Size() const { return Channels.size(); }
    const std::map<TMS_ChannelId, Entry>& GetChannels() const { return Channels; }
    static std::string ToString(Status status);

  private:
    std::map<TMS_ChannelId, Entry> Channels;
};

#endif
