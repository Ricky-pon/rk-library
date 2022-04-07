#ifndef RK_GRAPH_TEMPLATE_HPP
#define RK_GRAPH_TEMPLATE_HPP

#include <vector>

namespace rklib {

template <class T = int>
struct Edge {
    int from, to;
    T cost;
    int idx;

    Edge() = default;

    Edge(int from, int to, T cost = 1, int idx = 0)
        : from(from), to(to), cost(cost), idx(idx) {}
};

template <class T = int>
struct Graph {
    std::vector<std::vector<Edge<T>>> es;
    int edge_num;

    Graph(int n) : edge_num(0) { es.resize(n); }

    int size() { return es.size(); }

    void add_edge(int from, int to, T cost, int idx) {
        es[from].emplace_back(from, to, cost, idx);
        ++edge_num;
    }

    std::vector<Edge<T>> edges() {
        std::vector<Edge<T>> res(edge_num);
        for (int v = 0; v < (int)es.size(); ++v) {
            for (auto& e : es[v]) {
                res[e.idx] = e;
            }
        }
        return res;
    }

    std::vector<Edge<T>>& operator[](int i) { return es[i]; }
};

}  // namespace rklib

#endif  // RK_GRAPH_TEMPLATE_HPP