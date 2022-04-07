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
T div_floor(T a, T b) {
    if (b < 0) a *= -1, b *= -1;
    return a >= 0 ? a / b : (a + 1) / b - 1;
}

template <class T>
T div_ceil(T a, T b) {
    if (b < 0) a *= -1, b *= -1;
    return a > 0 ? (a - 1) / b + 1 : a / b;
}

template <typename T>
struct CoordComp {
    std::vector<T> v;
    bool sorted;

    CoordComp() : sorted(false) {}

    int size() { return v.size(); }

    void add(T x) { v.push_back(x); }

    void build() {
        std::sort(v.begin(), v.end());

        v.erase(std::unique(v.begin(), v.end()), v.end());
        sorted = true;
    }

    int get_idx(T x) {
        assert(sorted);
        return lower_bound(v.begin(), v.end(), x) - v.begin();
    }

    T &operator[](int i) { return v[i]; }
};

}  // namespace rklib

#endif  // RK_UTILITY_HPP