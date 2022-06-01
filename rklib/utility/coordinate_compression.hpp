#ifndef RK_COORD_COMP_HPP
#define RK_COORD_COMP_HPP

namespace rklib {

template <typename T = int>
struct CoordComp {
    std::vector<T> v;

    CoordComp() {}
    CoordComp(std::vector<T> &a) : v(a) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    int size() { return v.size(); }

    int get_idx(T x) { return lower_bound(v.begin(), v.end(), x) - v.begin(); }

    T &operator[](int i) { return v[i]; }
};

}  // namespace rklib

#endif  // RK_COORD_COMP_HPP