#ifndef RK_KRUSCAL_HPP
#define RK_KRUSCAL_HPP

#include <atcoder/dsu>
#include <rklib/graph/graph_template.hpp>
#include <vector>

namespace rklib {

template <class T>
std::pair<T, std::vector<int>> kruscal(Graph<T> &gr) {
    auto es = gr.edges();
    std::sort(es.begin(), es.end(), [](const Edge<T> &a, const Edge<T> &b) {
        return a.cost < b.cost;
    });
    atcoder::dsu uf(gr.size());
    T cost = 0;
    std::vector<int> mst_es;
    for (auto &e : es) {
        if (uf.same(e.from, e.to)) continue;
        uf.merge(e.from, e.to);
        cost += e.cost;
        mst_es.push_back(e.idx);
    }
    return {cost, mst_es};
}

}  // namespace rklib

#endif  // RK_KRUSCAL_HPP