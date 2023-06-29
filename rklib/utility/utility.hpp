#ifndef RK_UTILITY_HPP
#define RK_UTILITY_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace rklib {

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
    os << "[";
    for (auto d : v) os << d << ", ";
    return os << "]";
}

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

}  // namespace rklib

#endif  // RK_UTILITY_HPP