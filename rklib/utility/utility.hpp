#ifndef RK_UTILITY_HPP
#define RK_UTILITY_HPP

#include <algorithm>
#include <cassert>
#include <vector>

namespace rklib {

template <class T>
bool chmax(T &a, const T &b) {
    if (a < b) {
        a = b;
        return true;
    }
    return false;
}

template <class T>
bool chmin(T &a, const T &b) {
    if (a > b) {
        a = b;
        return true;
    }
    return false;
}

template <class T>
bool chmin_non_negative(T &a, const T &b) {
    if (a < 0 || a > b) {
        a = b;
        return true;
    }
    return false;
}

template <class T>
T div_floor(T a, T b) {
    if (b < 0) a *= -1, b *= -1;
    return a >= 0 ? a / b : (a + 1) / b - 1;
}

template <class T>
T div_ceil(T a, T b) {
    if (b < 0) a *= -1, b *= -1;
    return a > 0 ? (a - 1) / b + 1 : a / b;
}

}  // namespace rklib

#endif  // RK_UTILITY_HPP