#ifndef RK_MIN_OF_MOD_OF_LINEAR_HPP
#define RK_MIN_OF_MOD_OF_LINEAR_HPP

#include <cassert>
#include <numeric>
#include <rklib/utility/utility.hpp>
#include <vector>

namespace rklib {

// (x_0, x_1, ..., x_n)
// (d_0, d_1, ..., d_{n-1})
std::pair<std::vector<long long>, std::vector<long long>> min_mod_segments(
    long long m, long long a, long long b) {
    assert(0 <= a && a < m);
    assert(0 <= b && b < m);

    std::vector<long long> xs = {0};
    std::vector<long long> ds;
    auto g = std::gcd(a, m);
    a /= g, b /= g, m /= g;

    long long p = 0, q = 1, r = 1, s = 1;
    long long det_l = m - a, det_r = a;
    long long x = 0, y = b;

    while (y) {
        int k = det_r / det_l;
        det_r %= det_l;
        if (det_r == 0) {
            --k;
            det_r = det_l;
        }
        r += k * p;
        s += k * q;
        while (1) {
            int k = std::max(0LL, div_ceil(det_l - y, det_r));
            if (det_l - k * det_r <= 0) break;
            det_l -= k * det_r;
            p += k * r;
            q += k * s;

            k = y / det_l;
            y -= k * det_l;
            x += q * k;
            xs.push_back(x);
            ds.push_back(q);
        }
        k = det_l / det_r;
        det_l -= k * det_r;
        p += k * r;
        q += k * s;
        assert(std::min({p, q, r, s}) >= 0);
    }
    return {xs, ds};
}

// (a * x + b mod m, x)
std::pair<long long, long long> min_mod(long long n, long long m, long long a,
                                        long long b) {
    auto [xs, ds] = min_mod_segments(m, a, b);
    for (size_t i = 0; i < ds.size(); ++i) {
        if (xs[i + 1] <= n - 1) continue;
        auto k = (n - 1 - xs[i]) / ds[i];
        auto x = xs[i] + ds[i] * k;
        return {(a * x + b) % m, x};
    }
    return {(a * xs.back() + b) % m, xs.back()};
}

}  // namespace rklib

#endif  // RK_MIN_OF_MOD_OF_LINEAR_HPP