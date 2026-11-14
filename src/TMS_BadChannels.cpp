#include "TMS_BadChannels.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "toml.hpp"

#include "TMS_Manager.h"

std::string TMS_BadChannels::ToString(Status status) {
  switch (status) {
    case Status::Missing: return "missing";
    case Status::Dead: return "dead";
    default: return "good";
  }
}

TMS_BadChannels::TMS_BadChannels(const std::string& filename) {
  if (filename.empty()) return;

  std::string path = filename;
  if (path[0] != '/') {
    if (!std::getenv("TMS_DIR")) throw std::runtime_error("TMS_DIR must be set to find the bad channel file " + filename);
    path = std::string(std::getenv("TMS_DIR")) + "/" + filename;
  }

  const auto data = toml::parse(path);
  if (!data.contains("BadChannel")) return;
  for (const auto& table : toml::find<std::vector<toml::table>>(data, "BadChannel")) {
    const toml::value entry(table);
    const TMS_ChannelId channel{toml::find<int>(entry, "Plane"), toml::find<int>(entry, "Bar"),
                                toml::find<int>(entry, "View")};
    const std::string state = toml::find<std::string>(entry, "Status");
    Entry e;
    if (state == "missing") e.State = Status::Missing;
    else if (state == "dead") e.State = Status::Dead;
    else throw std::runtime_error("Bad channel file " + path + ": Status must be \"missing\" or \"dead\", not \"" + state + "\"");
    e.Reason = toml::find_or<std::string>(entry, "Reason", "");
    if (!Channels.emplace(channel, e).second)
      throw std::runtime_error("Bad channel file " + path + ": channel listed twice (plane " + std::to_string(channel.Plane) +
                               ", bar " + std::to_string(channel.Bar) + ", view " + std::to_string(channel.View) + ")");
  }
}

TMS_BadChannels::Status TMS_BadChannels::GetStatus(const TMS_ChannelId& channel) const {
  const auto it = Channels.find(channel);
  return it == Channels.end() ? Status::Good : it->second.State;
}

const TMS_BadChannels& TMS_BadChannels::GetInstance() {
  static const TMS_BadChannels instance = [] {
    const std::string file = TMS_Manager::GetInstance().Get_Decoder_BadChannelFile();
    TMS_BadChannels list(file);
    std::cout << "Bad channel list \"" << file << "\": " << list.Size() << " channels" << std::endl;
    return list;
  }();
  return instance;
}
