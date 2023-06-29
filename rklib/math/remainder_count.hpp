#ifndef RK_REMAINDER_COUNT_HPP
#define RK_REMAINDER_COUNT_HPP

#include <algorithm>

namespace rklib {

namespace internal {

template <class T>
T remainder_count(T r, T b, T m) {
    return r / m * b + std::min(b, r % m);
}

}  // namespace internal

// Number of integer x s.t.
// - x in [l, r)
// - x mod m in [a, b)
template <class T>
T remainder_count(T l, T r, T a, T b, T m) {
    assert(l >= 0);
    assert(m >= 1);

    if (l >= r || a >= b) return 0;
    if (m <= a || b < 0) return 0;
    chmax(a, T(0));
    chmin(b, m);

    auto res = internal::remainder_count(r, b, m);
    if (l >= 1) res -= internal::remainder_count(l, b, m);
    if (a >= 1) res -= internal::remainder_count(r, a, m);
    if (l >= 1 && a >= 1) res += internal::remainder_count(l, a, m);

    return res;
}

}  // namespace rklib

#endif  // RK_REMAINDER_COUNT_HPP