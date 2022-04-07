#ifndef RK_TREE_DIAMETER_HPP
#define RK_TREE_DIAMETER_HPP

#include <rklib/graph/graph_template.hpp>
#include <rklib/utility/utility.hpp>

namespace rklib {

template <class T>
std::tuple<T, std::vector<int>, std::vector<int>> tree_diameter(Graph<T> &gr) {
    T weight = 0;
    std::vector<int> path_v, path_e;
    auto dfs = [&](auto self, int v, int pv, int pe, T w) -> bool {
        bool flag = false;
        if (chmax(weight, w)) {
            flag = true;
            path_v.clear();
            path_e.clear();
        }
        for (auto &e : gr[v]) {
            if (e.to == pv) continue;
            flag |= self(self, e.to, v, e.idx, w + e.cost);
        }
        if (flag) {
            path_v.push_back(v);
            if (pe >= 0) path_e.push_back(pe);
        }
        return flag;
    };
    dfs(dfs, 0, -1, -1, 0);
    int s = path_v.front();
    weight = 0;
    path_v.clear();
    path_e.clear();
    dfs(dfs, s, -1, -1, 0);
    std::reverse(path_v.begin(), path_v.end());
    std::reverse(path_e.begin(), path_e.end());
    return {weight, path_v, path_e};
}

}  // namespace rklib

#endif  // RK_TREE_DIAMETER_HPP