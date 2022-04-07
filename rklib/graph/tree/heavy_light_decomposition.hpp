#ifndef RK_HEAVY_LIGHT_DECOMPOSITION_HPP
#define RK_HEAVY_LIGHT_DECOMPOSITION_HPP

#include <rklib/graph/graph_template.hpp>
#include <rklib/utility/utility.hpp>

namespace rklib {

struct HeavyLightDecomposition {
   public:
    std::vector<int> idx, par, head, out;

    HeavyLightDecomposition(Graph<> gr, int root = 0) {
        idx.resize(gr.size());
        par.resize(gr.size());
        head.resize(gr.size());
        out.resize(gr.size());

        sz_dfs(gr, root, -1);
        int cur = 0;
        idx[root] = cur++;
        head[root] = root;
        hld_dfs(gr, root, -1, root, cur);
    }

    int lca(int u, int v) {
        while (true) {
            if (idx[u] > idx[v]) std::swap(u, v);
            if (head[u] == head[v]) return u;
            v = par[head[v]];
        }
    }

    std::pair<std::vector<std::pair<int, int>>,
              std::vector<std::pair<int, int>>>
    get_path(int u, int v) {
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
    int sz_dfs(Graph<> &gr, int v, int pv) {
        int ret = 1, max_sz = 0;
        for (int i = 0; i < (int)gr[v].size(); ++i) {
            if (gr[v][i].to == pv) continue;
            int tmp = sz_dfs(gr, gr[v][i].to, v);
            ret += tmp;
            if (rklib::chmax(max_sz, tmp)) std::swap(gr[v][0], gr[v][i]);
        }
        return ret;
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