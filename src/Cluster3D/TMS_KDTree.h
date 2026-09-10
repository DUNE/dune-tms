#ifndef __TMS_KDTREE_H__
#define __TMS_KDTREE_H__

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

// Minimal static 3D KD-tree: build once from a fixed point array, then answer
// radius queries. Not a general-purpose spatial index -- no insertion,
// deletion, or k-NN support -- just what TMS_SpacePointDBScan needs for its
// epsilon-neighborhood lookups.
struct TMS_KDNode {
  int point_index = -1;
  int left = -1;
  int right = -1;
};

class TMS_KDTree {
  public:
    // pts: physical (x,y,z) positions, one per space point, indexed the same
    // way as the caller's space-point array. The tree stores indices into
    // this array, not copies of the coordinates, so `pts` must outlive the
    // tree.
    explicit TMS_KDTree(const std::vector<std::array<double, 3>> &pts) : _pts(pts) {
      if (_pts.empty()) return;
      std::vector<int> indices(_pts.size());
      std::iota(indices.begin(), indices.end(), 0);
      _pool.reserve(_pts.size());
      _root = Build(indices, 0, static_cast<int>(indices.size()), 0);
    }

    // Indices (into the array passed to the constructor) of all points within
    // `radius` of `query_point`, INCLUDING an exact match if present.
    std::vector<int> RadiusQuery(const std::array<double, 3> &query_point, double radius) const {
      std::vector<int> out;
      RadiusSearch(_root, query_point, radius, out, 0);
      return out;
    }

    // Convenience overload: query around an existing point in the tree by its index.
    std::vector<int> RadiusQuery(int query_index, double radius) const {
      return RadiusQuery(_pts[query_index], radius);
    }

  private:
    // Recursively median-split indices[lo,hi) on axis (depth % 3), building
    // nodes into the flat _pool. Returns the pool index of the subtree root,
    // or -1 for an empty range.
    int Build(std::vector<int> &indices, int lo, int hi, int depth) {
      if (lo >= hi) return -1;

      int axis = depth % 3;
      int mid = lo + (hi - lo) / 2;
      std::nth_element(indices.begin() + lo, indices.begin() + mid, indices.begin() + hi,
                        [&](int a, int b) { return _pts[a][axis] < _pts[b][axis]; });

      TMS_KDNode node;
      node.point_index = indices[mid];
      int node_pool_index = static_cast<int>(_pool.size());
      _pool.push_back(node);

      int left = Build(indices, lo, mid, depth + 1);
      int right = Build(indices, mid + 1, hi, depth + 1);
      _pool[node_pool_index].left = left;
      _pool[node_pool_index].right = right;

      return node_pool_index;
    }

    // Standard KD-tree radius search: test this node's point, then descend
    // into the near child (the side the query point lies on) always, and into
    // the far child only if the query's distance to the splitting plane is
    // within radius (i.e. the far side could still contain points within
    // range). `depth` re-derives the split axis, matching Build()'s convention.
    void RadiusSearch(int node_index, const std::array<double, 3> &q, double radius,
                       std::vector<int> &out, int depth) const {
      if (node_index < 0) return;
      const TMS_KDNode &node = _pool[node_index];
      const std::array<double, 3> &p = _pts[node.point_index];

      double dx = p[0] - q[0];
      double dy = p[1] - q[1];
      double dz = p[2] - q[2];
      double dist2 = dx * dx + dy * dy + dz * dz;
      if (dist2 <= radius * radius) out.push_back(node.point_index);

      int axis = depth % 3;
      double diff = q[axis] - p[axis];

      int near_child = (diff < 0.0) ? node.left : node.right;
      int far_child = (diff < 0.0) ? node.right : node.left;

      RadiusSearch(near_child, q, radius, out, depth + 1);
      if (std::fabs(diff) <= radius) {
        RadiusSearch(far_child, q, radius, out, depth + 1);
      }
    }

    const std::vector<std::array<double, 3>> &_pts;
    std::vector<TMS_KDNode> _pool;
    int _root = -1;
};

#endif
