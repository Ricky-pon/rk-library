#ifndef RK_DIJKSTRA_HPP
#define RK_DIJKSTRA_HPP

#include <limits>
#include <queue>
#include <rklib/graph/graph_template.hpp>
#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

template <class T>
std::pair<std::vector<T>, std::vector<int>> dijkstra(Graph<T> &gr, int s) {
    constexpr auto INF = std::numeric_limits<T>::max();
    std::vector<T> dist(gr.size(), INF);
    std::vector<int> pre(gr.size(), -1);
    std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>,
                        std::greater<std::pair<T, int>>>
        que;
    que.emplace(0, s);
    dist[s] = 0;
    while (!que.empty()) {
        auto [d, v] = que.top();
        que.pop();
        if (d > dist[v]) continue;
        for (auto &e : gr[v])
            if (chmin(dist[e.to], dist[v] + e.cost)) {
                que.emplace(dist[e.to], e.to);
                pre[e.to] = e.idx;
            }
    }
    return {dist, pre};
}

template <class T>
std::pair<T, std::vector<int>> dijkstra(Graph<T> &gr, int s, int t) {
    auto [dist, pre] = dijkstra(gr, s);
    std::vector<int> res;
    auto es = gr.edges();
    for (int e = pre[t]; e >= 0; e = pre[es[e].from]) {
        res.push_back(e);
    }
    reverse(res.begin(), res.end());
    return {dist[t], res};
}

}  // namespace rklib

#endif  // RK_DIJKSTRA_HPP