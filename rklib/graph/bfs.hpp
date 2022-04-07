#ifndef RK_BFS_HPP
#define RK_BFS_HPP

#include <limits>
#include <queue>
#include <rklib/graph/graph_template.hpp>
#include <vector>

namespace rklib {

template <class T>
std::pair<std::vector<T>, std::vector<int>> bfs(Graph<T> &gr, int s) {
    constexpr auto INF = std::numeric_limits<T>::max();
    std::vector<T> dist(gr.size(), INF);
    std::vector<int> pre(gr.size(), -1);
    std::queue<int> que;
    que.push(s);
    dist[s] = 0;
    while (!que.empty()) {
        auto v = que.front();
        que.pop();
        for (auto &e : gr[v])
            if (dist[e.to] == INF) {
                dist[e.to] = dist[v] + e.cost;
                que.push(e.to);
                pre[e.to] = e.idx;
            }
    }
    return {dist, pre};
}

template <class T>
std::pair<T, std::vector<int>> bfs(Graph<T> &gr, int s, int t) {
    auto [dist, pre] = bfs(gr, s);
    std::vector<int> res;
    auto es = gr.edges();
    for (int e = pre[t]; e >= 0; e = pre[es[e].from]) {
        res.push_back(e);
    }
    reverse(res.begin(), res.end());
    return {dist[t], res};
}

}  // namespace rklib

#endif  // RK_BFS_HPP