#ifndef RK_UTILITY_HPP
#define RK_UTILITY_HPP

#include <algorithm>
#include <cassert>
#include <vector>

namespace rklib {

// a <- max(a, b)
template <class T>
bool chmax(T &a, const T &b) {
    if (a < b) {
        a = b;
        return true;
    }
    return false;
}

// a <- min(a, b)
template <class T>
bool chmin(T &a, const T &b) {
    if (a > b) {
        a = b;
        return true;
    }
    return false;
}

// if a < 0: a <- b
// else: a <- min(a, b)
template <class T>
bool chmin_non_negative(T &a, const T &b) {
    if (a < 0 || a > b) {
        a = b;
        return true;
    }
    return false;
}

// floor(num / den)
template <class T>
T div_floor(T num, T den) {
    if (den < 0) num = -num, den = -den;
    return num >= 0 ? num / den : (num + 1) / den - 1;
}

// ceil(num / den)
template <class T>
T div_ceil(T num, T den) {
    if (den < 0) num = -num, den = -den;
    return num <= 0 ? num / den : (num - 1) / den + 1;
}

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

#endif  // RK_UTILITY_HPP