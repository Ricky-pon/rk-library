#ifndef RK_HEAVY_LIGHT_DECOMPOSITION_HPP
#define RK_HEAVY_LIGHT_DECOMPOSITION_HPP

#include <cassert>
#include <rklib/graph/graph_template.hpp>
#include <rklib/utility/utility.hpp>

namespace rklib {

struct HeavyLightDecomposition {
   public:
    int n;
    std::vector<int> idx, par, head, out, dep;

    HeavyLightDecomposition(Graph<> gr, int root = 0) : n(gr.size()) {
        idx.resize(gr.size());
        par.resize(gr.size());
        head.resize(gr.size());
        out.resize(gr.size());
        dep.resize(gr.size());

        size_dfs(gr, root, -1, 0);
        int cur = 0;
        idx[root] = cur++;
        head[root] = root;
        hld_dfs(gr, root, -1, root, cur);
    }

    int position(int v) {
        assert(0 <= v && v < n);
        return idx[v];
    }

    std::vector<int> node_sequence() {
        std::vector<int> res(idx.size());
        for (size_t v = 0; v < idx.size(); ++v) {
            res[idx[v]] = v;
        }
        return res;
    }

    std::pair<int, int> subtree_range(int v) {
        assert(0 <= v && v < n);
        return {idx[v], out[v]};
    }

    int lca(int u, int v) {
        assert(0 <= u && u < n);
        assert(0 <= v && v < n);
        while (true) {
            if (idx[u] > idx[v]) std::swap(u, v);
            if (head[u] == head[v]) return u;
            v = par[head[v]];
        }
    }

    int distance(int u, int v) { return dep[u] + dep[v] - dep[lca(u, v)] * 2; }

    int kth_ancestor_position(int v, int k) {
        assert(0 <= v && v < n);
        assert(0 <= k);
        while (true) {
            if (int d = idx[v] - idx[head[v]] + 1; d <= k) {
                k -= d;
                v = par[head[v]];
            } else {
                return idx[v] - k;
            }
        }
    }

    int jump_on_tree_position(int s, int t, int d) {
        assert(0 <= s && s < n);
        assert(0 <= t && t < n);
        assert(0 <= d);
        int v = lca(s, t), dist_st = dep[s] + dep[t] - dep[v] * 2;
        if (dist_st < d) return -1;
        if (dep[s] - dep[v] >= d) {
            return kth_ancestor_position(s, d);
        }
        return kth_ancestor_position(t, dist_st - d);
    }

    std::pair<std::vector<std::pair<int, int>>,
              std::vector<std::pair<int, int>>>
    path_ranges(int u, int v) {
        assert(0 <= u && u < n);
        assert(0 <= v && v < n);
        std::vector<std::pair<int, int>> up, down;
        while (true) {
            if (head[u] == head[v]) {
                if (idx[u] < idx[v])
                    down.emplace_back(idx[u], idx[v] + 1);
                else
                    up.emplace_back(idx[v], idx[u] + 1);
                std::reverse(up.begin(), up.end());
                std::reverse(down.begin(), down.end());
                return {up, down};
            } else {
                if (idx[u] < idx[v]) {
                    down.emplace_back(idx[head[v]], idx[v] + 1);
                    v = par[head[v]];
                } else {
                    up.emplace_back(idx[head[u]], idx[u] + 1);
                    u = par[head[u]];
                }
            }
        }
    }

   private:
    int size_dfs(Graph<> &gr, int v, int pv, int d) {
        int res = 1, max_sz = 0;
        dep[v] = d;
        for (int i = 0; i < (int)gr[v].size(); ++i) {
            if (gr[v][i].to == pv) continue;
            int tmp = size_dfs(gr, gr[v][i].to, v, d + 1);
            res += tmp;
            if (rklib::chmax(max_sz, tmp)) std::swap(gr[v][0], gr[v][i]);
        }
        return res;
    }

    void hld_dfs(Graph<> &gr, int v, int pv, int h, int &cur) {
        for (auto &e : gr[v]) {
            if (e.to == pv) continue;
            idx[e.to] = cur++;
            par[e.to] = v;
            head[e.to] = (gr[v][0].to == e.to ? h : e.to);
            hld_dfs(gr, e.to, v, head[e.to], cur);
        }
        out[v] = cur;
    }
};

}  // namespace rklib

#endif  // RK_HEAVY_LIGHT_DECOMPOSITION_HPP