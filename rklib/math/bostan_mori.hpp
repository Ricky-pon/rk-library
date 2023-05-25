#ifndef RK_BOSTAN_MORI_HPP
#define RK_BOSTAN_MORI_HPP

#include <atcoder/convolution>
#include <vector>

namespace rklib {

template <class T>
T bostan_mori_fps(std::vector<T> p, std::vector<T> q, long long n) {
    while (n >= 1) {
        auto q_minus = q;
        for (size_t i = 1; i < q.size(); i += 2) q_minus[i] = -q[i];
        auto pq = atcoder::convolution(p, q_minus);
        auto v = atcoder::convolution(q, q_minus);
        p.clear();
        q.clear();
        for (size_t i = 0; i < v.size(); i += 2) q.push_back(v[i]);
        if (n % 2 == 0) {
            for (size_t i = 0; i < pq.size(); i += 2) p.push_back(pq[i]);
        } else {
            for (size_t i = 1; i < pq.size(); i += 2) p.push_back(pq[i]);
        }
        n >>= 1;
    }
    return p[0];
}

template <class T>
T bostan_mori_rec(std::vector<T> a, std::vector<T> c, long long n) {
    int d = a.size();
    if (n < (long long)d) {
        return a[n];
    }
    std::vector<T> q(d + 1);
    q[0] = 1;
    for (int i = 0; i < d; i++) {
        q[i + 1] = -c[i];
    }
    a = atcoder::convolution(a, q);
    a.resize(d);
    return bostan_mori_fps(a, q, n);
}

}  // namespace rklib

#endif  // RK_BOSTAN_MORI_HPP