#pragma once

#include <stdexcept>

constexpr int TMS_VertexIdScale = 1000000;

inline long long TMS_MakeGlobalVertexID(int run_id, int vertex_id) {
  if (vertex_id >= TMS_VertexIdScale) {
    throw std::runtime_error("Vertex ID exceeds the global vertex ID encoding range");
  }
  return static_cast<long long>(run_id) * TMS_VertexIdScale + static_cast<long long>(vertex_id);
}
